import test from 'node:test';
import assert from 'node:assert/strict';
import {
 defaults3d,init3d,step3d,drive3d,tracerMap3d,
 triangleProximity,periodicProximity,surfaceCandidates,cluster3d,diagnostics3d
} from '../strips3d.mjs';

const finite3=v=>v.every(Number.isFinite);
const approximately=(a,b,eps=1e-7)=>Math.abs(a-b)<eps;

test('default 3D mesh has actual finite-width triangulated surfaces',()=>{
 const s=init3d({strips:8});
 assert.equal(s.strips.length,8);
 assert.equal(s.vertices.length,8*defaults3d.cols*defaults3d.rows);
 assert.equal(s.triangles.length,8*2*(defaults3d.cols-1)*(defaults3d.rows-1));
 assert.ok(s.vertices.every(v=>finite3(v.p)&&finite3(v.v)));
});

test('same seed yields deterministic 3D initialization',()=>{
 const a=init3d({strips:5,seed:71}),b=init3d({strips:5,seed:71});
 assert.deepEqual(a.vertices.map(v=>v.p),b.vertices.map(v=>v.p));
});

test('a triangle pair has nonzero separation when parallel and zero when intersecting',()=>{
 const A=[[0,0,0],[1,0,0],[0,1,0]];
 const B=A.map(p=>[p[0],p[1],.12]);
 assert.ok(approximately(triangleProximity(A,B).distance,.12));
 const C=[[.25,.25,-1],[.25,.25,1],[.8,.75,0]];
 assert.ok(triangleProximity(A,C).distance<1e-9);
});

test('periodic triangle proximity detects surfaces on opposite sides of box',()=>{
 const s={p:{box:7},vertices:[],triangles:[]};
 const P=[[.05,1,1],[.45,1,1],[.05,1.5,1]];
 const Q=[[6.95,1,1.07],[6.55,1,1.07],[6.95,1.5,1.07]];
 for(const x of P.concat(Q))s.vertices.push({p:x});
 const a={ids:[0,1,2]},b={ids:[3,4,5]};
 assert.ok(periodicProximity(s,a,b).distance<.16);
});

test('broad phase includes a forced inter-strip overlap',()=>{
 const s=init3d({strips:2});
 const a=s.strips[0].ids,b=s.strips[1].ids;
 for(let i=0;i<a.length;i++)s.vertices[b[i]].p=s.vertices[a[i]].p.map((x,k)=>x+(k===2?.1:0));
 assert.ok(surfaceCandidates(s).some(([i,j])=>s.triangles[i].strip!==s.triangles[j].strip));
 step3d(s,1);
 assert.equal(s.links.size,1);
 assert.ok(s.memory.size>0);
 assert.ok(s.maxOverlap>0);
 const d=diagnostics3d(s);
 assert.equal(d.fraction,1);
 assert.equal(d.edges,1);
});

test('3D tangential contact memory is cleared when the surfaces separate',()=>{
 const s=init3d({strips:2});
 const a=s.strips[0].ids,b=s.strips[1].ids;
 for(let i=0;i<a.length;i++){
  s.vertices[b[i]].p=s.vertices[a[i]].p.map((x,k)=>x+(k===2?.1:0));
  s.vertices[b[i]].v=[.15,0,0];
 }
 step3d(s,2);
 assert.ok(s.memory.size>0);
 const size=s.memory.size;
 for(const id of b)s.vertices[id].p=s.vertices[id].p.map(x=>x+3.35);
 step3d(s,1);
 assert.equal(s.memory.size,0);
 assert.equal(s.links.size,0);
 assert.ok(s.completed.length>0);
});

test('three consecutive shears have an exact volume-preserving passive tracer map',()=>{
 const P={box:7,K:.38,tau:2.5};
 const initial=[1.1,2.4,3.2],x=initial.slice();
 for(let phase=0;phase<3;phase++){const u=drive3d(x,(phase+.1)*P.tau,P);for(let k=0;k<3;k++)x[k]+=P.tau*u[k];}
 const map=tracerMap3d(initial,P.K,P.box);
 for(let k=0;k<3;k++)assert.ok(approximately(x[k],map[k]));
 // The Jacobian is an upper/lower shear composition, determinant 1.
 const eps=1e-5,columns=[];
 for(let j=0;j<3;j++){
  const plus=initial.slice(),minus=initial.slice();plus[j]+=eps;minus[j]-=eps;
  const a=tracerMap3d(plus,P.K,P.box),b=tracerMap3d(minus,P.K,P.box);
  columns.push(a.map((v,k)=>(v-b[k])/(2*eps)));
 }
 const [u,v,w]=columns;
 const det=u[0]*(v[1]*w[2]-w[1]*v[2])-v[0]*(u[1]*w[2]-w[1]*u[2])+w[0]*(u[1]*v[2]-v[1]*u[2]);
 assert.ok(Math.abs(det-1)<1e-5);
});

test('three-dimensional driven elastic dynamics remain finite for soft and stiff strips',()=>{
 for(const bend of [0,14,60]){
  const s=init3d({strips:6,bending:bend,seed:19});
  step3d(s,120);
  assert.ok(s.vertices.every(v=>finite3(v.p)&&finite3(v.v)&&finite3(v.f)));
  const d=diagnostics3d(s);
  assert.ok(d.fraction>=1/6&&d.fraction<=1);
  assert.ok(Number.isFinite(d.maxOverlap));
 }
});

test('the largest contact cluster is empty for no strips, and isolated for disconnected strips',()=>{
 const a=init3d({strips:0}),b=init3d({strips:4});
 assert.equal(cluster3d(a).fraction,0);
 assert.equal(cluster3d(b).largest,1);
 assert.equal(cluster3d(b).fraction,.25);
});
