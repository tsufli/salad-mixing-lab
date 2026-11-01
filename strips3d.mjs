// Exploratory 3D periodic deformable-strip DEM; no lentils, adhesion or fluid feedback.
// Each strip is a triangulated elastic spring-grid surface. Contact uses closest
// triangle features (vertex-face and edge-edge), with persistent tangential springs.
export const defaults3d={
 strips:22,length:3.2,width:.9,cols:7,rows:3,box:8,bending:14,
 friction:.75,K:.5,tau:2.5,dt:.002,seed:17,
 spring:280,normalStiffness:650,normalDamp:10,tangentialStiffness:110,
 tangentialDamp:1.5,drag:7,thickness:.16
};
const pi=Math.PI;
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const mod=(x,L)=>((x%L)+L)%L;
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]];
const scale=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>Math.sqrt(dot(a,a));
export function seededRandom(seed){let s=seed>>>0;return()=>{s=(1664525*s+1013904223)>>>0;return s/4294967296;};}
export function drive3d(pos,t,p){
 const {box:L,K,tau}=p,q=L*K/tau,phase=Math.floor((t+1e-10)/tau)%3;
 if(phase===0)return[q*Math.sin(2*pi*pos[1]/L),0,0];
 if(phase===1)return[0,q*Math.sin(2*pi*pos[2]/L),0];
 return[0,0,q*Math.sin(2*pi*pos[0]/L)];
}
export function tracerMap3d(r,K,L=8){
 let x=mod(r[0]/L+K*Math.sin(2*pi*r[1]/L),1);
 let y=mod(r[1]/L+K*Math.sin(2*pi*r[2]/L),1);
 let z=mod(r[2]/L+K*Math.sin(2*pi*x),1);
 return [x*L,y*L,z*L];
}
function randomFrame(rand){
 const u=rand(),v=rand(),w=rand(),a=Math.sqrt(1-u),b=Math.sqrt(u);
 const q=[a*Math.sin(2*pi*v),a*Math.cos(2*pi*v),b*Math.sin(2*pi*w),b*Math.cos(2*pi*w)];
 const [x,y,z,r]=q;
 return [[1-2*(y*y+z*z),2*(x*y+z*r),2*(x*z-y*r)],
 [2*(x*y-z*r),1-2*(x*x+z*z),2*(y*z+x*r)]];
}
function makeStrip(s,index,center,frame){
 const {p,vertices,triangles,springs,bends,strips}=s,{cols:C,rows:R,length,width}=p,ids=[];
 for(let i=0;i<C;i++){const row=[];for(let j=0;j<R;j++){
  const a=(i/(C-1)-.5)*length,b=(j/(R-1)-.5)*width;
  const x=center.map((v,k)=>v+a*frame[0][k]+b*frame[1][k]);
  row.push(vertices.length);vertices.push({p:x,v:[0,0,0],f:[0,0,0],m:.22,strip:index});
 }ids.push(row);}
 function spring(a,b,k=1){const d=norm(sub(vertices[a].p,vertices[b].p));springs.push({a,b,rest:d,k});}
 function bend(indices,coef,k=1){bends.push({indices,coef,k});}
 for(let i=0;i<C;i++)for(let j=0;j<R;j++){
  if(i+1<C)spring(ids[i][j],ids[i+1][j],1);
  if(j+1<R)spring(ids[i][j],ids[i][j+1],1);
  if(i+1<C&&j+1<R){
   spring(ids[i][j],ids[i+1][j+1],.65);
   spring(ids[i+1][j],ids[i][j+1],.65);
   bend([ids[i][j],ids[i+1][j],ids[i][j+1],ids[i+1][j+1]],[1,-1,-1,1],.4);
  }
  if(i>0&&i+1<C)bend([ids[i-1][j],ids[i][j],ids[i+1][j]],[1,-2,1],1);
  if(j>0&&j+1<R)bend([ids[i][j-1],ids[i][j],ids[i][j+1]],[1,-2,1],.45);
  if(i+1<C&&j+1<R){
   const a=ids[i][j],b=ids[i+1][j],c=ids[i][j+1],d=ids[i+1][j+1];
   triangles.push({id:triangles.length,strip:index,ids:[a,b,c],i,j});
   triangles.push({id:triangles.length,strip:index,ids:[b,d,c],i,j});
  }
 }
 strips.push({id:index,ids:ids.flat(),grid:ids,color:index});
}
function triPoints(s,t){return t.ids.map(i=>s.vertices[i].p);}
function center3(p){return scale(add(add(p[0],p[1]),p[2]),1/3);}
function radius3(p,c){return Math.max(...p.map(x=>norm(sub(x,c))));}
function closestTri(p,a,b,c){
 const ab=sub(b,a),ac=sub(c,a),ap=sub(p,a);
 let d1=dot(ab,ap),d2=dot(ac,ap);
 if(d1<=0&&d2<=0)return{p:a,w:[1,0,0]};
 const bp=sub(p,b),d3=dot(ab,bp),d4=dot(ac,bp);
 if(d3>=0&&d4<=d3)return{p:b,w:[0,1,0]};
 let vc=d1*d4-d3*d2;
 if(vc<=0&&d1>=0&&d3<=0){let v=d1/(d1-d3);return{p:add(a,scale(ab,v)),w:[1-v,v,0]};}
 const cp=sub(p,c),d5=dot(ab,cp),d6=dot(ac,cp);
 if(d6>=0&&d5<=d6)return{p:c,w:[0,0,1]};
 let vb=d5*d2-d1*d6;
 if(vb<=0&&d2>=0&&d6<=0){let w=d2/(d2-d6);return{p:add(a,scale(ac,w)),w:[1-w,0,w]};}
 let va=d3*d6-d5*d4;
 if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){let w=(d4-d3)/((d4-d3)+(d5-d6));return{p:add(b,scale(sub(c,b),w)),w:[0,1-w,w]};}
 const denom=va+vb+vc;
 if(Math.abs(denom)<1e-16)return{p:a,w:[1,0,0]};
 let v=vb/denom,w=vc/denom;return{p:add(a,add(scale(ab,v),scale(ac,w))),w:[1-v-w,v,w]};
}
function closestSegments(a,b,c,d){
 const u=sub(b,a),v=sub(d,c),w=sub(a,c),A=dot(u,u),B=dot(u,v),C=dot(v,v),D=dot(u,w),E=dot(v,w),den=A*C-B*B;
 let s=0,t=0;
 if(A<1e-14&&C<1e-14)return{s,t};
 if(A<1e-14){t=clamp(E/C,0,1);return{s,t};}
 if(C<1e-14){s=clamp(-D/A,0,1);return{s,t};}
 if(den>1e-14)s=clamp((B*E-C*D)/den,0,1);
 t=clamp((B*s+E)/C,0,1);
 s=clamp((B*t-D)/A,0,1);
 t=clamp((B*s+E)/C,0,1);
 return{s,t};
}
// Closest surface features: 6 vertex/face and 9 edge/edge tests.
function segmentTriangle(P,Q,T){
 const dir=sub(Q,P),e1=sub(T[1],T[0]),e2=sub(T[2],T[0]),h=cross(dir,e2),det=dot(e1,h);
 if(Math.abs(det)<1e-12)return null;
 const inv=1/det,sp=sub(P,T[0]),u=inv*dot(sp,h);
 if(u<0||u>1)return null;
 const q=cross(sp,e1),v=inv*dot(dir,q);
 if(v<0||u+v>1)return null;
 const t=inv*dot(e2,q);
 if(t<0||t>1)return null;
 return{t,bary:[1-u-v,u,v],point:add(P,scale(dir,t))};
}
export function triangleProximity(A,B){
 let best=Infinity,wa=null,wb=null,pa=null,pb=null;
 function offer(x,y,a,b){const d2=dot(sub(x,y),sub(x,y));if(d2<best){best=d2;wa=a;wb=b;pa=x;pb=y;}}
 for(let i=0;i<3;i++){const a=[0,0,0];a[i]=1;const q=closestTri(A[i],...B);offer(A[i],q.p,a,q.w);}
 for(let j=0;j<3;j++){const b=[0,0,0];b[j]=1;const q=closestTri(B[j],...A);offer(q.p,B[j],q.w,b);}
 const E=[[0,1],[1,2],[2,0]];
 for(const [i,j] of E)for(const [k,l] of E){
  const {s,t}=closestSegments(A[i],A[j],B[k],B[l]);
  const x=add(A[i],scale(sub(A[j],A[i]),s)),y=add(B[k],scale(sub(B[l],B[k]),t));
  const a=[0,0,0],b=[0,0,0];a[i]=1-s;a[j]=s;b[k]=1-t;b[l]=t;offer(x,y,a,b);
 }
 // Edge through face interior is an intersection, even when neither edge
 // endpoint is close to the opposite triangle (a gap in naive proximity tests).
 for(const [i,j] of E){
  const q=segmentTriangle(A[i],A[j],B);
  if(q){const a=[0,0,0];a[i]=1-q.t;a[j]=q.t;offer(q.point,q.point,a,q.bary);}
  const r=segmentTriangle(B[i],B[j],A);
  if(r){const b=[0,0,0];b[i]=1-r.t;b[j]=r.t;offer(r.point,r.point,r.bary,b);}
 }
 return{distance:Math.sqrt(best),wa,wb,pa,pb};
}
export function periodicProximity(s,a,b){
 const A=triPoints(s,a),B=triPoints(s,b),ca=center3(A),cb=center3(B),L=s.p.box;
 const shift=ca.map((x,k)=>Math.round((x-cb[k])/L)*L);
 const shifted=B.map(v=>add(v,shift));
 return triangleProximity(A,shifted);
}
export function init3d(options={}){
 const p={...defaults3d,...options},rng=seededRandom(p.seed),L=p.box;
 const s={p,t:0,vertices:[],triangles:[],springs:[],bends:[],strips:[],memory:new Map(),pairAges:new Map(),completed:[],links:new Set(),maxOverlap:0,contacts:0,placementFailed:0};
 for(let i=0;i<p.strips;i++){
  let placed=false;
  for(let trial=0;trial<130;trial++){
   const center=[rng()*L,rng()*L,rng()*L],frame=randomFrame(rng);
   const nv=s.vertices.length,nt=s.triangles.length,ns=s.springs.length,nb=s.bends.length;
   makeStrip(s,i,center,frame);
   let overlap=false;
   const fresh=s.triangles.slice(nt);
   for(const a of fresh){if(overlap)break;const P=triPoints(s,a),ca=center3(P),ra=radius3(P,ca);
    for(let j=0;j<nt;j++){const b=s.triangles[j],Q=triPoints(s,b),cb=center3(Q),rb=radius3(Q,cb);
     const dist=Math.hypot(...ca.map((x,k)=>x-cb[k]-Math.round((x-cb[k])/L)*L));
     if(dist>ra+rb+p.thickness*1.15)continue;
     if(periodicProximity(s,a,b).distance<p.thickness*1.15){overlap=true;break;}
    }
   }
   if(!overlap){placed=true;break;}
   s.vertices.length=nv;s.triangles.length=nt;s.springs.length=ns;s.bends.length=nb;s.strips.length=i;
  }
  if(!placed){s.placementFailed=p.strips-i;break;}
 }
 return s;
}
function collisionCandidates(s){
 const {triangles,p}=s,L=p.box,info=new Array(triangles.length);
 let maxRadius=0;
 for(let i=0;i<triangles.length;i++){
  const P=triPoints(s,triangles[i]),c=center3(P),r=radius3(P,c);
  maxRadius=Math.max(maxRadius,r);
  info[i]={c,r,lo:[0,1,2].map(k=>Math.min(...P.map(v=>v[k]))),
    hi:[0,1,2].map(k=>Math.max(...P.map(v=>v[k])))};
 }
 // Cell width >= sum of any two triangle radii plus contact thickness.
 // Thus every true contact is found in this cell or its 26 neighbors,
 // even if the deformable triangles stretch beyond their nominal size.
 const n=Math.max(2,Math.floor(L/(2*maxRadius+p.thickness+1e-8))),cell=L/n,hash=new Map();
 const key=(a,b,c)=>(a*n+b)*n+c;
 for(let i=0;i<triangles.length;i++){
  const q=info[i].c.map(x=>Math.floor(mod(x,L)/cell)%n);info[i].q=q;
  const id=key(...q);if(!hash.has(id))hash.set(id,[]);hash.get(id).push(i);
 }
 const out=[];
 for(let i=0;i<triangles.length;i++){
  const a=triangles[i],v=info[i],seen=new Set();
  for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++){
   const cellIndex=key(mod(v.q[0]+dx,n),mod(v.q[1]+dy,n),mod(v.q[2]+dz,n));
   for(const j of hash.get(cellIndex)||[]){
    if(j<=i||seen.has(j))continue;seen.add(j);
    const b=triangles[j];
    if(a.strip===b.strip&&Math.abs(a.i-b.i)<=1)continue;
    const w=info[j],dist=Math.hypot(...v.c.map((x,k)=>x-w.c[k]-Math.round((x-w.c[k])/L)*L));
    if(dist>v.r+w.r+p.thickness)continue;
    const shift=v.c.map((x,k)=>Math.round((x-w.c[k])/L)*L);
    let nearby=true;
    for(let k=0;k<3;k++)if(v.lo[k]-(w.hi[k]+shift[k])>p.thickness||(w.lo[k]+shift[k])-v.hi[k]>p.thickness){nearby=false;break;}
    if(nearby)out.push([i,j]);
   }
  }
 }
 return out;
}
export function surfaceCandidates(s){return collisionCandidates(s);}
function applyElastic(s){
 const {vertices:v,p}=s;
 for(const e of s.springs){
  const a=v[e.a],b=v[e.b],d=sub(b.p,a.p),length=Math.max(1e-12,norm(d)),
  force=scale(d,p.spring*e.k*(length-e.rest)/length);
  for(let k=0;k<3;k++){a.f[k]+=force[k];b.f[k]-=force[k];}
 }
 for(const e of s.bends){
  const x=[0,0,0];for(let i=0;i<e.indices.length;i++){const q=v[e.indices[i]].p;for(let k=0;k<3;k++)x[k]+=e.coef[i]*q[k];}
  for(let i=0;i<e.indices.length;i++){const a=v[e.indices[i]],c=-p.bending*e.k*e.coef[i];for(let k=0;k<3;k++)a.f[k]+=c*x[k];}
 }
}
function contactForce(s,ta,tb,id,dt){
 const {p,vertices:v}=s,prox=periodicProximity(s,ta,tb);
 if(prox.distance>=p.thickness)return false;
 let normal;
 if(prox.distance>1e-8)normal=scale(sub(prox.pa,prox.pb),1/prox.distance);
 else{
  // The normal is ambiguous for intersecting surfaces. Choose the face normal
  // with sign oriented by their centers; this is only a soft-penalty fallback.
  const A=triPoints(s,ta),B=triPoints(s,tb),n=cross(sub(A[1],A[0]),sub(A[2],A[0]));
  const m=norm(n);normal=m>1e-12?scale(n,1/m):[1,0,0];
  if(dot(normal,sub(center3(A),center3(B)))<0)normal=scale(normal,-1);
 }
 const wa=prox.wa,wb=prox.wb,va=[0,0,0],vb=[0,0,0];
 for(let z=0;z<3;z++)for(let k=0;k<3;k++){
  va[k]+=wa[z]*v[ta.ids[z]].v[k];vb[k]+=wb[z]*v[tb.ids[z]].v[k];
 }
 const relative=sub(va,vb),vn=dot(relative,normal),vt=sub(relative,scale(normal,vn));
 const fn=Math.max(0,p.normalStiffness*(p.thickness-prox.distance)-p.normalDamp*vn);
 const old=s.memory.get(id);
 const memory=old?sub(old.xi,scale(normal,dot(old.xi,normal))):[0,0,0];
 let xi=add(memory,scale(vt,dt));
 let Ft=add(scale(xi,-p.tangentialStiffness),scale(vt,-p.tangentialDamp));
 const mag=norm(Ft),limit=p.friction*fn;
 if(mag>limit){Ft=mag>0?scale(Ft,limit/mag):[0,0,0];xi=p.tangentialStiffness?scale(add(Ft,scale(vt,p.tangentialDamp)),-1/p.tangentialStiffness):[0,0,0];}
 const F=add(scale(normal,fn),Ft);
 for(let z=0;z<3;z++)for(let k=0;k<3;k++){
  v[ta.ids[z]].f[k]+=wa[z]*F[k];v[tb.ids[z]].f[k]-=wb[z]*F[k];
 }
 s.memory.set(id,{xi,age:(old?.age||0)+dt,stripA:ta.strip,stripB:tb.strip});
 s.maxOverlap=Math.max(s.maxOverlap,p.thickness-prox.distance);
 s.contacts++;
 if(ta.strip!==tb.strip){
  const a=Math.min(ta.strip,tb.strip),b=Math.max(ta.strip,tb.strip);
  s.links.add(a+':'+b);
 }
 return true;
}
function recenterStrips(s){
 const L=s.p.box;
 for(const strip of s.strips){
  const mean=[0,0,0];
  for(const id of strip.ids)for(let k=0;k<3;k++)mean[k]+=s.vertices[id].p[k]/strip.ids.length;
  const shifts=mean.map(x=>Math.floor(x/L)*L);
  if(shifts.some(x=>x!==0))for(const id of strip.ids)for(let k=0;k<3;k++)s.vertices[id].p[k]-=shifts[k];
 }
}
export function step3d(s,n=1){
 const{p,vertices}=s;
 for(let step=0;step<n;step++){
  for(const a of vertices){const u=drive3d(a.p,s.t,p);for(let k=0;k<3;k++)a.f[k]=p.drag*a.m*(u[k]-a.v[k]);}
  applyElastic(s);
  s.contacts=0;s.maxOverlap=0;s.links=new Set();
  const active=new Set();
  for(const [i,j] of collisionCandidates(s)){
   const key=i+':'+j;if(contactForce(s,s.triangles[i],s.triangles[j],key,p.dt))active.add(key);
  }
  for(const [key,item] of s.memory)if(!active.has(key))s.memory.delete(key);
  for(const [key,age] of s.pairAges)if(!s.links.has(key)){s.completed.push(age);s.pairAges.delete(key);}
  for(const key of s.links)s.pairAges.set(key,(s.pairAges.get(key)||0)+p.dt);
  if(s.completed.length>400)s.completed.splice(0,s.completed.length-400);
  for(const a of vertices)for(let k=0;k<3;k++){
   a.v[k]=clamp(a.v[k]+p.dt*a.f[k]/a.m,-8,8);
   a.p[k]+=p.dt*a.v[k];
  }
  recenterStrips(s);s.t+=p.dt;
 }
 return s;
}
export function cluster3d(s){
 const N=s.strips.length;
 if(!N)return{largest:0,fraction:0,degree:0,meanAge:0,edges:0,component:[]};
 const parents=Array.from({length:N},(_,i)=>i),adj=Array.from({length:N},()=>new Set());
 const root=i=>{while(parents[i]!==i){parents[i]=parents[parents[i]];i=parents[i];}return i;};
 for(const edge of s.links){const[a,b]=edge.split(':').map(Number);parents[root(a)]=root(b);adj[a].add(b);adj[b].add(a);}
 const groups=new Map();for(let i=0;i<N;i++){const r=root(i);if(!groups.has(r))groups.set(r,[]);groups.get(r).push(i);}
 let component=[];for(const g of groups.values())if(g.length>component.length)component=g;
 let sum=0;for(const value of s.pairAges.values())sum+=value;
 return{largest:component.length,fraction:component.length/N,degree:s.links.size*2/N,
 meanAge:s.pairAges.size?sum/s.pairAges.size:0,edges:s.links.size,component};
}
function solve3(A,b){
 const m=A.map((row,i)=>[...row,b[i]]);
 for(let i=0;i<3;i++){let pivot=i;for(let j=i+1;j<3;j++)if(Math.abs(m[j][i])>Math.abs(m[pivot][i]))pivot=j;
  if(Math.abs(m[pivot][i])<1e-10)return null;
  [m[i],m[pivot]]=[m[pivot],m[i]];
  const div=m[i][i];for(let k=i;k<4;k++)m[i][k]/=div;
  for(let j=0;j<3;j++)if(j!==i){const f=m[j][i];for(let k=i;k<4;k++)m[j][k]-=f*m[i][k];}
 }
 return m.map(r=>r[3]);
}
// Fraction of COM velocity variance explained by a single angular velocity,
// after subtracting common translation; near unity = rigid-body-like motion.
export function coherence3d(s,component=cluster3d(s).component){
 if(component.length<3)return null;
 const L=s.p.box,positions=[],velocities=[];
 for(const i of component){
  const strip=s.strips[i],P=[0,0,0],V=[0,0,0];
  for(const id of strip.ids)for(let k=0;k<3;k++){P[k]+=s.vertices[id].p[k]/strip.ids.length;V[k]+=s.vertices[id].v[k]/strip.ids.length;}
  positions.push(P);velocities.push(V);
 }
 const anchor=positions[0];
 const relative=positions.map(p=>p.map((x,k)=>anchor[k]+(x-anchor[k]-Math.round((x-anchor[k])/L)*L)));
 const meanP=[0,0,0],meanV=[0,0,0];
 for(let i=0;i<relative.length;i++)for(let k=0;k<3;k++){meanP[k]+=relative[i][k]/relative.length;meanV[k]+=velocities[i][k]/relative.length;}
 const A=[[0,0,0],[0,0,0],[0,0,0]],b=[0,0,0];let variance=0;
 for(let i=0;i<relative.length;i++){
  const r=sub(relative[i],meanP),v=sub(velocities[i],meanV),rr=dot(r,r),rv=cross(r,v);variance+=dot(v,v);
  for(let j=0;j<3;j++){b[j]+=rv[j];for(let k=0;k<3;k++)A[j][k]+=(j===k?rr:0)-r[j]*r[k];}
 }
 if(variance<1e-9)return 1;
 const omega=solve3(A,b);if(!omega)return null;
 let residual=0;for(let i=0;i<relative.length;i++){const r=sub(relative[i],meanP),v=sub(velocities[i],meanV),d=sub(v,cross(omega,r));residual+=dot(d,d);}
 return clamp(1-residual/variance,0,1);
}
export function diagnostics3d(s){
 const c=cluster3d(s),avgFinished=s.completed.length?s.completed.reduce((a,b)=>a+b,0)/s.completed.length:0;
 return{t:s.t,fraction:c.fraction,edges:c.edges,degree:c.degree,meanAge:c.meanAge,
 completedAge:avgFinished,coherence:coherence3d(s,c.component),maxOverlap:s.maxOverlap,contacts:s.contacts,strips:s.strips.length};
}
