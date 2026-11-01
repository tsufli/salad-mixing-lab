// 2D periodic, flexible-fiber and disk DEM. Dimensionless exploratory model.
export const defaults={fibers:30,disks:180,segments:5,aspect:6,bending:8,friction:.6,K:.7,initial:'uniform',halfPeriod:2.5,dt:.002,kn:650,etaN:10,kt:130,etaT:2,ks:400,drag:8,seed:17};
const wrap=(x,L)=>((x%L)+L)%L;
const clamp=(v,a,b)=>Math.min(b,Math.max(a,v));
export function rng(seed){let s=seed>>>0;return()=>{s=(1664525*s+1013904223)>>>0;return s/4294967296;};}
export function drivingVelocity(x,y,t,p,W=18,H=12){const tau=p.halfPeriod||2.5,phase=Math.floor((t+1e-9)/tau)%2;return phase===0?{x:W*p.K/tau*Math.sin(2*Math.PI*y/H),y:0}:{x:0,y:H*p.K/tau*Math.sin(2*Math.PI*x/W)};}
export function tracerMap(x,y,K,W=18,H=12){const X=wrap(x/W+K*Math.sin(2*Math.PI*y/H),1),Y=wrap(y/H+K*Math.sin(2*Math.PI*X),1);return{x:W*X,y:H*Y};}
function delta(x,L){return x-Math.round(x/L)*L;}
export function init(options={}){
 const p={...defaults,...options},rand=rng(p.seed),W=18,H=12,particles=[],fibers=[];
 function safe(x,y,r,margin=.04){for(const a of particles){if(Math.hypot(delta(x-a.x,W),delta(y-a.y,H))<r+a.r+margin)return false;}return true;}
 for(let f=0;f<p.fibers;f++){
  const length=.23*p.aspect,l=length/(p.segments-1),r=.105;
  let chosen=null;
  for(let a=0;a<1800;a++){
   const x=rand()*W,y=rand()*H,theta=rand()*Math.PI*2,ux=Math.cos(theta),uy=Math.sin(theta);
   const nodes=Array.from({length:p.segments},(_,i)=>({x:wrap(x+(i-(p.segments-1)/2)*l*ux,W),y:wrap(y+(i-(p.segments-1)/2)*l*uy,H)}));
   if(nodes.every(n=>safe(n.x,n.y,r,.01))){chosen=nodes;break;}
  }
  if(!chosen)break;
  const ids=[];
  for(const n of chosen){ids.push(particles.length);particles.push({type:'node',fiber:f,x:n.x,y:n.y,ux:n.x,uy:n.y,vx:0,vy:0,fx:0,fy:0,m:.38,r});}
  fibers.push({ids,l});
 }
 for(let i=0;i<p.disks;i++){
  for(let j=0;j<2500;j++){
   const x=rand()*W,y=rand()*H,r=.09;
   if(p.initial==='left'&&x>W/2)continue;
   if(!safe(x,y,r)){continue;}
   particles.push({type:'disk',fiber:-1,x,y,ux:x,uy:y,vx:0,vy:0,fx:0,fy:0,m:.25,r});break;
  }
 }
 return{p,W,H,particles,fibers,t:0,contacts:0,rodContacts:[],memory:new Map(),ended:[],contactClock:0,observations:0};
}
export function candidatePairs(s){
 const{particles,W,H}=s,nx=18,ny=12,cells=Array.from({length:nx*ny},()=>[]);
 for(let i=0;i<particles.length;i++){const a=particles[i],cx=Math.floor(a.x/W*nx)%nx,cy=Math.floor(a.y/H*ny)%ny;cells[cy*nx+cx].push(i);}
 const pairs=[];
 for(let i=0;i<particles.length;i++){
  const a=particles[i],cx=Math.floor(a.x/W*nx)%nx,cy=Math.floor(a.y/H*ny)%ny;
  for(let dy=-1;dy<=1;dy++)for(let dx=-1;dx<=1;dx++){
   for(const j of cells[wrap(cy+dy,ny)*nx+wrap(cx+dx,nx)])if(j>i&&!(a.fiber>=0&&a.fiber===particles[j].fiber))pairs.push([i,j]);
  }
 }
 return pairs;
}
function elasticForces(s){
 const{particles,p,W,H}=s;
 for(const fiber of s.fibers){
  const{ids,l}=fiber;
  for(let i=0;i<ids.length-1;i++){
   const a=particles[ids[i]],b=particles[ids[i+1]];
   const dx=delta(b.x-a.x,W),dy=delta(b.y-a.y,H),d=Math.max(1e-10,Math.hypot(dx,dy));
   const force=p.ks*(d-l),fx=force*dx/d,fy=force*dy/d;
   a.fx+=fx;a.fy+=fy;b.fx-=fx;b.fy-=fy;
  }
  // Discrete bending energy E_b=(kappa/2)sum_i|r_{i-1}-2r_i+r_{i+1}|^2.
  // Derivative is evaluated analytically; neighboring displacements are minimum-image.
  for(let i=1;i<ids.length-1;i++){
   const a=particles[ids[i-1]],b=particles[ids[i]],c=particles[ids[i+1]];
   const ex=delta(c.x-b.x,W)-delta(b.x-a.x,W);
   const ey=delta(c.y-b.y,H)-delta(b.y-a.y,H);
   a.fx-=p.bending*ex;a.fy-=p.bending*ey;
   b.fx+=2*p.bending*ex;b.fy+=2*p.bending*ey;
   c.fx-=p.bending*ex;c.fy-=p.bending*ey;
  }
 }
}
export function contact(s,i,j,dt=s.p.dt){
 const{particles,p,W,H}=s,a=particles[i],b=particles[j];
 const dx=delta(b.x-a.x,W),dy=delta(b.y-a.y,H),d=Math.hypot(dx,dy);
 const key=i+':'+j,old=s.memory.get(key);
 if(d>=a.r+b.r){return null;}
 const nx=d>1e-10?dx/d:1,ny=d>1e-10?dy/d:0,tx=-ny,ty=nx;
 const dvx=b.vx-a.vx,dvy=b.vy-a.vy,vn=dvx*nx+dvy*ny,vt=dvx*tx+dvy*ty;
 const fn=Math.max(0,p.kn*(a.r+b.r-d)-p.etaN*vn);
 // Tangential history is a scalar for 2D contacts; old contact displacement is
 // projected onto the new tangent if the normal rotates.
 let xi=(old?old.xi*(old.tx*tx+old.ty*ty):0)+vt*dt;
 let ft=-p.kt*xi-p.etaT*vt;
 const limit=p.friction*fn;
 if(Math.abs(ft)>limit){
   ft=clamp(ft,-limit,limit);
   xi=p.kt>0?(-ft-p.etaT*vt)/p.kt:0;
 }
 const fx=fn*nx+ft*tx,fy=fn*ny+ft*ty;
 a.fx-=fx;a.fy-=fy;b.fx+=fx;b.fy+=fy;
 const age=(old?old.age:0)+dt;
 s.memory.set(key,{xi,tx,ty,age,fi:a.fiber,fj:b.fiber});
 if(a.fiber>=0&&b.fiber>=0&&a.fiber!==b.fiber)s.rodContacts.push([a.fiber,b.fiber,age]);
 return key;
}
export function step(s,n=1){
 const{particles,p,W,H}=s;
 for(let z=0;z<n;z++){
  for(const a of particles){
   const v=drivingVelocity(a.x,a.y,s.t,p,W,H);
   a.fx=p.drag*a.m*(v.x-a.vx);
   a.fy=p.drag*a.m*(v.y-a.vy);
  }
  elasticForces(s);
  const active=new Set();s.rodContacts=[];s.contacts=0;
  for(const [i,j] of candidatePairs(s)){const key=contact(s,i,j,p.dt);if(key!==null){active.add(key);s.contacts++;}}
  for(const [key,item] of s.memory)if(!active.has(key)){s.ended.push(item.age);s.memory.delete(key);}
  if(s.ended.length>600)s.ended.splice(0,s.ended.length-600);
  for(const a of particles){
   a.vx=clamp(a.vx+p.dt*a.fx/a.m,-12,12);a.vy=clamp(a.vy+p.dt*a.fy/a.m,-12,12);
   a.ux+=p.dt*a.vx;a.uy+=p.dt*a.vy;
   a.x=wrap(a.ux,W);a.y=wrap(a.uy,H);
  }
  s.t+=p.dt;
 }
 return s;
}
export function clusterStats(s){
 const N=s.fibers.length;if(!N)return{largest:0,largestFraction:0,meanDegree:0,edges:0,meanAge:0};
 const parent=Array.from({length:N},(_,i)=>i),neighbors=Array.from({length:N},()=>new Set());
 const root=i=>{while(parent[i]!==i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
 let age=0;
 for(const [a,b,t] of s.rodContacts){parent[root(a)]=root(b);neighbors[a].add(b);neighbors[b].add(a);age+=t;}
 const sizes=new Map();for(let i=0;i<N;i++){const a=root(i);sizes.set(a,(sizes.get(a)||0)+1);}
 const edges=neighbors.reduce((sum,v)=>sum+v.size,0)/2;
 return{largest:Math.max(...sizes.values()),largestFraction:Math.max(...sizes.values())/N,meanDegree:2*edges/N,edges,meanAge:s.rodContacts.length?age/s.rodContacts.length:0};
}
export function structureFactor(s,maxMode=5){
 const disks=s.particles.filter(a=>a.type==='disk'),N=disks.length;
 if(!N)return[];
 const groups=new Map(),dk=.25;
 for(let nx=-maxMode;nx<=maxMode;nx++)for(let ny=-maxMode;ny<=maxMode;ny++){
  if(!nx&&!ny)continue;
  const kx=2*Math.PI*nx/s.W,ky=2*Math.PI*ny/s.H,k=Math.hypot(kx,ky),bucket=Math.round(k/dk);
  let re=0,im=0;for(const a of disks){let phase=kx*a.x+ky*a.y;re+=Math.cos(phase);im+=Math.sin(phase);}
  const v=groups.get(bucket)||{k:0,S:0,n:0};v.k+=k;v.S+=(re*re+im*im)/N;v.n++;groups.set(bucket,v);
 }
 return[...groups.values()].map(v=>({k:v.k/v.n,S:v.S/v.n})).sort((a,b)=>a.k-b.k);
}
export function mixing(s){
 const counts=new Array(12).fill(0),disks=s.particles.filter(a=>a.type==='disk');
 for(const a of disks)counts[clamp(Math.floor(a.x/s.W*4),0,3)+4*clamp(Math.floor(a.y/s.H*3),0,2)]++;
 const mean=disks.length/12;return counts.reduce((v,c)=>v+(c-mean)**2,0)/(12*(mean*mean||1));
}
export function contactSurvival(s,threshold=1){
 const active=[...s.memory.values()].filter(x=>x.fi>=0&&x.fj>=0&&x.fi!==x.fj);
 return{active:active.length,olderThanThreshold:active.filter(x=>x.age>=threshold).length,meanAge:active.length?active.reduce((v,x)=>v+x.age,0)/active.length:0,meanFinished:s.ended.length?s.ended.reduce((v,x)=>v+x,0)/s.ended.length:0};
}
