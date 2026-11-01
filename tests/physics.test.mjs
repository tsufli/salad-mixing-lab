import test from 'node:test';
import assert from 'node:assert/strict';
import {init,step,drivingVelocity,tracerMap,contactGeometry,pairForce,mixing} from '../physics.mjs';
test('deterministic initialization',()=>assert.deepEqual(init(),init()));
test('periodic disk contact across boundary',()=>assert.ok(contactGeometry({type:'disk',x:.05,y:3,r:.15},{type:'disk',x:11.95,y:3,r:.15}).overlap>.19));
test('rod-disk periodic contact exerts equal and opposite forces',()=>{
 const a={type:'rod',x:11.8,y:4,theta:0,len:1,r:.1,vx:0,vy:0,omega:0,fx:0,fy:0,tau:0},b={type:'disk',x:.15,y:4.15,r:.1,vx:0,vy:0,omega:0,fx:0,fy:0,tau:0};
 assert.ok(contactGeometry(a,b).overlap>0);pairForce(a,b,{stiffness:900,damping:13,friction:.5});
 assert.ok(Math.abs(a.fx+b.fx)<1e-8);assert.ok(Math.abs(a.fy+b.fy)<1e-8);
});
test('exact passive tracer map matches two sine shear steps',()=>{
 const x=2,y=1,K=.3,p={K,halfPeriod:2.5},u=drivingVelocity(x,y,0,p),x1=x+u.x*2.5,v=drivingVelocity(x1,y,2.5+1e-6,p),y1=y+v.y*2.5,m=tracerMap(x,y,K);
 assert.ok(Math.abs(m.x-x1)<1e-9);assert.ok(Math.abs(m.y-y1)<1e-9);
});
test('tracer map preserves area locally',()=>{
 const x=2,y=1,k=.7,h=1e-5,f=tracerMap,x1=f(x+h,y,k),x0=f(x-h,y,k),y1=f(x,y+h,k),y0=f(x,y-h,k);
 const det=((x1.x-x0.x)*(y1.y-y0.y)-(x1.y-x0.y)*(y1.x-y0.x))/(4*h*h);
 assert.ok(Math.abs(det-1)<1e-6);
});
test('finite simulation with periodic position wrapping',()=>{
 const s=init({rods:8,disks:25});step(s,500);
 assert.ok(s.items.every(o=>[o.x,o.y,o.vx,o.vy,o.theta,o.omega].every(Number.isFinite)&&o.x>=0&&o.x<s.W&&o.y>=0&&o.y<s.H));
 assert.ok(Number.isFinite(mixing(s)));
});
