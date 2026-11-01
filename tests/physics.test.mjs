import test from 'node:test';
import assert from 'node:assert/strict';
import {init,step,contact,candidatePairs,clusterStats,structureFactor,contactSurvival,drivingVelocity,tracerMap} from '../physics.mjs';
test('deterministic initialization and counts',()=>{
 const a=init(),b=init();assert.deepEqual(a,b);
 assert.equal(a.fibers.length,30);
 assert.equal(a.particles.filter(x=>x.type==='disk').length,180);
 assert.equal(a.particles.filter(x=>x.type==='node').length,150);
});
test('periodic contacts conserve pairwise momentum',()=>{
 const s=init({fibers:0,disks:0}),a={type:'disk',fiber:-1,x:.05,y:3,r:.15,vx:0,vy:0,fx:0,fy:0,m:.25},b={...a,x:17.95,fx:0,fy:0};
 s.particles=[a,b];assert.ok(contact(s,0,1));assert.ok(Math.abs(a.fx+b.fx)<1e-10&&Math.abs(a.fy+b.fy)<1e-10);
});
test('tangential contact memory accumulates under sustained slip and clears on separation',()=>{
 const s=init({fibers:0,disks:0,friction:1}),a={type:'disk',fiber:-1,x:1,y:1,r:.15,vx:0,vy:0,fx:0,fy:0,m:.25},b={...a,x:1.25,vy:.15};
 s.particles=[a,b];
 contact(s,0,1,.01);const first=s.memory.get('0:1').xi;
 contact(s,0,1,.01);const second=s.memory.get('0:1').xi;
 assert.ok(second>first&&first>0);
 b.x=3;step(s,1);assert.equal(s.memory.size,0);
});
test('flexible spring and bending simulations stay finite',()=>{
 for(const bending of [0,8,45]){const s=init({fibers:12,disks:60,bending});step(s,1500);
 assert.ok(s.particles.every(a=>[a.x,a.y,a.vx,a.vy,a.ux,a.uy].every(Number.isFinite)));
 }
});
test('uniform and left-half initialization',()=>{
 const a=init({fibers:6,disks:80}),b=init({fibers:6,disks:80,initial:'left'});
 assert.ok(a.particles.some(x=>x.type==='disk'&&x.x>a.W/2));
 assert.ok(b.particles.filter(x=>x.type==='disk').every(x=>x.x<=b.W/2));
});
test('broad phase includes all overlapping particle pairs',()=>{
 const s=init({fibers:10,disks:80}),pairs=new Set(candidatePairs(s).map(([i,j])=>i+':'+j));
 for(let i=0;i<s.particles.length;i++)for(let j=i+1;j<s.particles.length;j++){
 const a=s.particles[i],b=s.particles[j];if(a.fiber>=0&&a.fiber===b.fiber)continue;
 const dx=b.x-a.x-Math.round((b.x-a.x)/s.W)*s.W,dy=b.y-a.y-Math.round((b.y-a.y)/s.H)*s.H;
 if(Math.hypot(dx,dy)<a.r+b.r)assert.ok(pairs.has(i+':'+j));
 }
});
test('contact clusters and Sdd statistics are finite',()=>{
 const s=init();step(s,250);const c=clusterStats(s),f=structureFactor(s),q=contactSurvival(s);
 assert.ok(c.largest>=1&&c.largest<=s.fibers.length);
 assert.ok(c.largestFraction<=1&&Number.isFinite(c.meanDegree));
 assert.ok(f.length>3&&f.every(x=>Number.isFinite(x.S)&&x.S>=0&&x.k>0));
 assert.ok(Number.isFinite(q.meanAge));
});
test('exact passive-tracer sine map has unit local Jacobian',()=>{
 const x=2,y=1,k=.7,h=1e-5,f=tracerMap;
 const a=f(x+h,y,k),b=f(x-h,y,k),c=f(x,y+h,k),d=f(x,y-h,k);
 const det=((a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x))/(4*h*h);
 assert.ok(Math.abs(det-1)<1e-6);
 const u=drivingVelocity(x,y,0,{K:k,halfPeriod:2.5});assert.ok(Number.isFinite(u.x));
});
