import {init3d,step3d,diagnostics3d} from './strips3d.mjs';

let sim=null,playing=false,batch=3,nextSample=0,token=0,stepRate=0;
function sendFrame(){
 if(!sim)return;
 const positions=new Float32Array(sim.vertices.length*3);
 for(let i=0;i<sim.vertices.length;i++)for(let k=0;k<3;k++)positions[i*3+k]=sim.vertices[i].p[k];
 const status=diagnostics3d(sim);status.stepsPerSecond=stepRate;
 postMessage({type:'frame',positions:positions.buffer,status},[positions.buffer]);
}
function sendSample(){
 const d=diagnostics3d(sim);
 postMessage({type:'sample',data:{t:d.t,fraction:d.fraction,coherence:d.coherence,
   edges:d.edges,meanAge:d.meanAge,completedAge:d.completedAge,maxOverlap:d.maxOverlap}});
}
function cycle(id){
 if(!playing||id!==token||!sim)return;
 try{
  const wasPreparing=sim.stage!=='stirring',start=Date.now();
  step3d(sim,batch);
  const elapsed=Math.max(1,Date.now()-start),instant=batch*1000/elapsed;
  stepRate=stepRate?stepRate*.8+instant*.2:instant;
  // Record only after compression AND relaxation have completed.
  if(sim.stage==='stirring'&&(wasPreparing||sim.t>=nextSample)){
   sendSample();nextSample=sim.t+.2;
  }
  sendFrame();
  setTimeout(()=>cycle(id),0);
 }catch(error){
  playing=false;postMessage({type:'error',message:String(error?.stack||error)});
 }
}
self.onmessage=e=>{
 const{cmd}=e.data;
 if(cmd==='reset'){
  token++;playing=false;
  try{
   sim=init3d(e.data.params);nextSample=0;stepRate=0;
   postMessage({type:'topology',box:sim.p.box,targetBox:sim.prepTargetBox,
    strips:sim.strips.map(x=>({id:x.id,ids:x.ids,grid:x.grid})),
    triangles:sim.triangles.map(x=>({strip:x.strip,ids:x.ids})),
    actual:sim.strips.length,requested:sim.p.strips,vertices:sim.vertices.length,faces:sim.triangles.length,
    status:diagnostics3d(sim)});
   sendFrame();
   if(sim.stage==='stirring'){sendSample();nextSample=.2;}
   playing=e.data.autoPlay!==false;
   if(playing)cycle(++token);
  }catch(error){postMessage({type:'error',message:String(error?.stack||error)});}
 }else if(cmd==='pause'){playing=false;token++;postMessage({type:'paused'});}
 else if(cmd==='play'&&sim&&!playing){playing=true;cycle(++token);}
 else if(cmd==='batch'){batch=Math.max(1,Math.min(12,Math.round(e.data.value||1)));}
 else if(cmd==='snapshot')sendFrame();
};