// Dimensionless periodic-box discrete-element model. All three contact species are active.
// The imposed flow is alternating sine shear, not a physical fluid solver.
export const defaults={rods:36,disks:220,aspect:6,initial:'uniform',friction:.55,K:.7,halfPeriod:2.5,dt:.003,stiffness:900,damping:13,drag:7,seed:17};
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const mod=(x,L)=>((x%L)+L)%L;
export function rng(seed){let x=seed>>>0;return ()=>{x=(1664525*x+1013904223)>>>0;return x/4294967296;};}
function closest(a,b,p){const dx=b.x-a.x,dy=b.y-a.y,t=clamp(((p.x-a.x)*dx+(p.y-a.y)*dy)/(dx*dx+dy*dy||1),0,1);return{x:a.x+t*dx,y:a.y+t*dy};}
function segmentClosest(a,b,c,d){const u={x:b.x-a.x,y:b.y-a.y},v={x:d.x-c.x,y:d.y-c.y},w={x:a.x-c.x,y:a.y-c.y};const A=u.x*u.x+u.y*u.y,B=u.x*v.x+u.y*v.y,C=v.x*v.x+v.y*v.y,D=u.x*w.x+u.y*w.y,E=v.x*w.x+v.y*w.y;const den=A*C-B*B;let s=den>1e-12?clamp((B*E-C*D)/den,0,1):0,t=clamp((B*s+E)/(C||1),0,1);s=clamp((B*t-D)/(A||1),0,1);t=clamp((B*s+E)/(C||1),0,1);return[{x:a.x+s*u.x,y:a.y+s*u.y},{x:c.x+t*v.x,y:c.y+t*v.y}];}
function endpoints(o){const dx=Math.cos(o.theta)*o.len/2,dy=Math.sin(o.theta)*o.len/2;return[{x:o.x-dx,y:o.y-dy},{x:o.x+dx,y:o.y+dy}];}
function geometricPair(a,b){
 let pa,pb;
 if(a.type==='disk'&&b.type==='disk'){pa={x:a.x,y:a.y};pb={x:b.x,y:b.y};}
 else if(a.type==='rod'&&b.type==='disk'){const[s,e]=endpoints(a);pa=closest(s,e,b);pb={x:b.x,y:b.y};}
 else if(a.type==='disk'&&b.type==='rod'){const[s,e]=endpoints(b);pa={x:a.x,y:a.y};pb=closest(s,e,a);}
 else{const[s,e]=endpoints(a),[u,v]=endpoints(b);[pa,pb]=segmentClosest(s,e,u,v);}
 const dx=pb.x-pa.x,dy=pb.y-pa.y,dist=Math.hypot(dx,dy);
 return{pa,pb,nx:dist>1e-9?dx/dist:1,ny:dist>1e-9?dy/dist:0,overlap:a.r+b.r-dist};
}
// Minimum-image periodic contact geometry. All modeled rods are shorter than
// half the shortest box dimension, so any overlapping pair uses this image.
export function contactGeometry(a,b,W=12,H=8){
 const bx=b.x+Math.round((a.x-b.x)/W)*W;
 const by=b.y+Math.round((a.y-b.y)/H)*H;
 const shifted={...b,x:bx,y:by};
 return {...geometricPair(a,shifted),bx,by};
}
export function pairForce(a,b,p,W=12,H=8){
 const c=contactGeometry(a,b,W,H);
 if(c.overlap<=0)return false;
 const ra={x:c.pa.x-a.x,y:c.pa.y-a.y},rb={x:c.pb.x-c.bx,y:c.pb.y-c.by};
 const va={x:a.vx-a.omega*ra.y,y:a.vy+a.omega*ra.x},vb={x:b.vx-b.omega*rb.y,y:b.vy+b.omega*rb.x};
 const dvx=vb.x-va.x,dvy=vb.y-va.y,vn=dvx*c.nx+dvy*c.ny;
 const fn=Math.max(0,p.stiffness*c.overlap-p.damping*vn),tx=-c.ny,ty=c.nx,vt=dvx*tx+dvy*ty;
 const ft=-p.friction*fn*Math.tanh(vt/.06),fx=fn*c.nx+ft*tx,fy=fn*c.ny+ft*ty;
 a.fx-=fx;a.fy-=fy;b.fx+=fx;b.fy+=fy;
 a.tau-=ra.x*fy-ra.y*fx;b.tau+=rb.x*fy-rb.y*fx;
 return true;
}
export function drivingVelocity(x,y,t,p,W=12,H=8){
 const phase=Math.floor((t+1e-10)/(p.halfPeriod||2.5))%2;
 if(phase===0)return{x:W*p.K/(p.halfPeriod||2.5)*Math.sin(2*Math.PI*y/H),y:0};
 return{x:0,y:H*p.K/(p.halfPeriod||2.5)*Math.sin(2*Math.PI*x/W)};
}
// Exact stroboscopic map for a passive tracer across a full period.
export function tracerMap(x,y,K,W=12,H=8){
 const X=mod(x/W+K*Math.sin(2*Math.PI*y/H),1);
 const Y=mod(y/H+K*Math.sin(2*Math.PI*X),1);
 return{x:W*X,y:H*Y};
}
export function init(options={}){
 const p={...defaults,...options},random=rng(p.seed),items=[],W=18,H=12;let id=0;
 function place(type,r,len){
  for(let tries=0;tries<2500;tries++){
   const theta=random()*Math.PI*2,x=random()*W,y=random()*H;
   if(type==='disk'&&p.initial==='left'&&x>W*.5)continue;
   const m=type==='rod'?2:.24;
   const obj={id:id++,type,x,y,theta,r,len,m,I:type==='rod'?m*(len*len/12+r*r/2):m*r*r/2,vx:0,vy:0,omega:0,fx:0,fy:0,tau:0};
   if(items.every(o=>contactGeometry(o,obj,W,H).overlap<-.015)){items.push(obj);return;}
   id--;
  }
 }
 for(let i=0;i<p.rods;i++)place('rod',.115,.23*p.aspect);
 for(let i=0;i<p.disks;i++)place('disk',.095,0);
 return{p,items,W,H,t:0,contacts:0};
}
// Center-based broad phase. The cell dimensions exceed the largest possible
// sum of two bounding radii at the supported aspect ratios (<= 9).
export function candidatePairs(s){
 const {items,W,H}=s,nx=Math.floor(W/2.5),ny=Math.floor(H/2.5),cells=Array.from({length:nx*ny},()=>[]);
 for(let i=0;i<items.length;i++){
  const o=items[i],cx=Math.floor(o.x/W*nx)%nx,cy=Math.floor(o.y/H*ny)%ny;
  cells[cy*nx+cx].push(i);
 }
 const pairs=[];
 for(let i=0;i<items.length;i++){
  const o=items[i],cx=Math.floor(o.x/W*nx)%nx,cy=Math.floor(o.y/H*ny)%ny;
  for(let dy=-1;dy<=1;dy++)for(let dx=-1;dx<=1;dx++){
   const ix=mod(cx+dx,nx),iy=mod(cy+dy,ny);
   for(const j of cells[iy*nx+ix])if(j>i)pairs.push([i,j]);
  }
 }
 return pairs;
}
// A contact graph based on instantaneous geometric overlap, distinct from
// accumulated contact history or persistent entanglement.
export function clusterStats(s){
 const rods=s.items.map((o,i)=>o.type==='rod'?i:-1).filter(i=>i>=0),N=rods.length;
 if(!N)return {largest:0,largestFraction:0,meanDegree:0,edges:0};
 const parent=Array.from({length:s.items.length},(_,i)=>i),degree=new Array(s.items.length).fill(0);
 function root(i){while(parent[i]!==i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
 let edges=0;
 for(const [i,j] of candidatePairs(s)){
  if(s.items[i].type!=='rod'||s.items[j].type!=='rod')continue;
  if(contactGeometry(s.items[i],s.items[j],s.W,s.H).overlap>0){
   parent[root(i)]=root(j);degree[i]++;degree[j]++;edges++;
  }
 }
 const counts=new Map();
 for(const i of rods){const r=root(i);counts.set(r,(counts.get(r)||0)+1);}
 const largest=Math.max(...counts.values());
 return {largest,largestFraction:largest/N,meanDegree:2*edges/N,edges};
}
// Radially averaged disk structure factor S_dd(k), excluding the q=0 mode.
// Snapshot values include the finite-N random baseline S~1.
export function structureFactor(s,maxMode=5){
 const disks=s.items.filter(o=>o.type==='disk'),N=disks.length;
 if(!N)return [];
 const buckets=new Map(),dk=.25;
 for(let nx=-maxMode;nx<=maxMode;nx++)for(let ny=-maxMode;ny<=maxMode;ny++){
  if(nx===0&&ny===0)continue;
  const kx=2*Math.PI*nx/s.W,ky=2*Math.PI*ny/s.H,k=Math.hypot(kx,ky);
  const shell=Math.round(k/dk);
  let re=0,im=0;
  for(const o of disks){const phi=kx*o.x+ky*o.y;re+=Math.cos(phi);im+=Math.sin(phi);}
  const b=buckets.get(shell)||{k:0,S:0,n:0};
  b.k+=k;b.S+=(re*re+im*im)/N;b.n++;buckets.set(shell,b);
 }
 return [...buckets.values()].map(b=>({k:b.k/b.n,S:b.S/b.n})).sort((a,b)=>a.k-b.k);
}
export function step(s,n=1){
 const{items,p,W,H}=s;
 for(let z=0;z<n;z++){
  s.contacts=0;
  for(const o of items){
   const v=drivingVelocity(o.x,o.y,s.t,p,W,H);
   o.fx=p.drag*o.m*(v.x-o.vx);o.fy=p.drag*o.m*(v.y-o.vy);
   let desiredOmega=0;
   if(o.type==='rod'){
    const ex=Math.cos(o.theta)*o.len/2,ey=Math.sin(o.theta)*o.len/2;
    const va=drivingVelocity(o.x-ex,o.y-ey,s.t,p,W,H),vb=drivingVelocity(o.x+ex,o.y+ey,s.t,p,W,H);
    desiredOmega=(-(vb.x-va.x)*Math.sin(o.theta)+(vb.y-va.y)*Math.cos(o.theta))/o.len;
   }
   o.tau=p.drag*o.I*(desiredOmega-o.omega);
  }
  for(const [i,j] of candidatePairs(s))if(pairForce(items[i],items[j],p,W,H))s.contacts++;
  for(const o of items){
   o.vx=clamp(o.vx+p.dt*o.fx/o.m,-20,20);o.vy=clamp(o.vy+p.dt*o.fy/o.m,-20,20);
   o.omega=clamp(o.omega+p.dt*o.tau/o.I,-25,25);
   o.x=mod(o.x+p.dt*o.vx,W);o.y=mod(o.y+p.dt*o.vy,H);o.theta+=p.dt*o.omega;
  }
  s.t+=p.dt;
 }
 return s;
}
export function mixing(s){
 const counts=new Array(12).fill(0),disks=s.items.filter(o=>o.type==='disk');
 for(const o of disks)counts[clamp(Math.floor(o.x/s.W*4),0,3)+4*clamp(Math.floor(o.y/s.H*3),0,2)]++;
 const avg=disks.length/12;
 return counts.reduce((v,c)=>v+(c-avg)**2,0)/12/(avg*avg||1);
}
