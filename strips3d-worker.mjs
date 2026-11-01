import {init3d,step3d,diagnostics3d} from './strips3d.mjs';
let sim=null,playing=false,batch=2,nextSample=.2,token=0;
function sendFrame(){
 if(!sim)return;
 const positions=new Float32Array(sim.vertices.length*3);
 for(let i=0;i<sim.vertices.length;i++)for(let k=0;k<3;k++)positions[i*3+k]=sim.vertices[i].p[k];
 const status=diagnostics3d(sim);
 postMessage({type:'frame',positions:positions.buffer,status},[positions.buffer]);
}
function sendSample(){
 const d=diagnostics3d(sim);
 postMessage({type:'sample',data:{t:d.t,fraction:d.fraction,coherence:d.coherence,
   edges:d.edges,meanAge:d.meanAge,completedAge:d.completedAge,maxOverlap:d.maxOverlap}});
}
function cycle(id){
 if(!playing||id!==token||!sim)return;
 step3d(sim,batch);
 if(sim.t>=nextSample){sendSample();nextSample=sim.t+.2;}
 sendFrame();
 setTimeout(()=>cycle(id),0);
}
self.onmessage=e=>{
 const {cmd}=e.data;
 if(cmd==='reset'){
  token++;playing=false;
  try{
   sim=init3d(e.data.params);
   nextSample=.2;
   postMessage({type:'topology',box:sim.p.box,strips:sim.strips.map(x=>({id:x.id,ids:x.ids,grid:x.grid})),
    triangles:sim.triangles.map(x=>({strip:x.strip,ids:x.ids})),
    actual:sim.strips.length,requested:sim.p.strips,vertices:sim.vertices.length,faces:sim.triangles.length,
    status:diagnostics3d(sim)});
   sendFrame();sendSample();
   playing=e.data.autoPlay!==false;
   if(playing)cycle(++token);
  }catch(err){postMessage({type:'error',message:String(err?.stack||err)});}
 }else if(cmd==='pause'){playing=false;token++;postMessage({type:'paused'});}
 else if(cmd==='play'&&sim&&!playing){playing=true;cycle(++token);}
 else if(cmd==='batch'){batch=Math.max(1,Math.min(5,Math.round(e.data.value||1)));}
 else if(cmd==='snapshot'){sendFrame();}
};