import { Canvas, useFrame } from "@react-three/fiber";
import { ContactShadows, Environment, Html, OrbitControls, RoundedBox, Sparkles } from "@react-three/drei";
import { Suspense, useMemo, useRef } from "react";
import * as THREE from "three";
import type { CompanyMarketData } from "../../features/market-war/types";
import { performance, territoryX, unitScale } from "../../features/market-war/calculations";
import { useMarketWarStore } from "../../stores/marketWarStore";

const BLUE="#3c8cff", RED="#ff594f", GOLD="#f6c866";
const unitPositions = (items:CompanyMarketData[]) => {
  const sides={USA:items.filter(c=>c.country==="USA"),CHINA:items.filter(c=>c.country==="CHINA")};
  const map=new Map<string,[number,number,number]>();
  for(const country of ["USA","CHINA"] as const) sides[country].forEach((c,i)=>{
    const col=i%5,row=Math.floor(i/5); const base=country==="USA"?-7.4:7.4;
    map.set(c.id,[base+(country==="USA"?1:-1)*(col%2)*1.3,.72,-4+col*2+row*.7]);
  }); return map;
};

function Tree({x,z,color}:{x:number;z:number;color:string}){return <group position={[x,.7,z]}><mesh><cylinderGeometry args={[.1,.14,.7,6]}/><meshStandardMaterial color="#70543b"/></mesh><mesh position={[0,.55,0]}><coneGeometry args={[.48,1,7]}/><meshStandardMaterial color={color}/></mesh></group>}
function WindTurbine({x,z}:{x:number;z:number}){return <group position={[x,.85,z]}><mesh><cylinderGeometry args={[.07,.12,1.6,8]}/><meshStandardMaterial color="#dce8ed"/></mesh><mesh position={[0,.8,0]} rotation={[0,0,.2]}><boxGeometry args={[.08,1.4,.08]}/><meshStandardMaterial color="#f5fbff"/></mesh><mesh position={[0,.8,0]} rotation={[0,0,2.3]}><boxGeometry args={[.08,1.4,.08]}/><meshStandardMaterial color="#f5fbff"/></mesh></group>}
function BattlefieldTerrain({frontline}:{frontline:number}){
  const trees=useMemo(()=>Array.from({length:26},(_,i)=>({x:(i%2?-1:1)*(3+(i*1.91)%8),z:-5.4+(i*2.67)%10})),[]);
  return <group>
    <mesh position={[0,-.45,0]}><boxGeometry args={[25,1,13]}/><meshStandardMaterial color="#17251f" roughness={.95}/></mesh>
    <mesh position={[-6,0.08,0]}><boxGeometry args={[12,.18,12]}/><meshStandardMaterial color="#285b72" roughness={.9}/></mesh>
    <mesh position={[6,0.08,0]}><boxGeometry args={[12,.18,12]}/><meshStandardMaterial color="#74413d" roughness={.9}/></mesh>
    <mesh position={[0,.18,0]} rotation={[-Math.PI/2,0,0]}><planeGeometry args={[1.2,12]}/><meshStandardMaterial color="#183e52" roughness={.25} metalness={.1}/></mesh>
    <mesh position={[-6,.2,0]} rotation={[-Math.PI/2,0,0]}><planeGeometry args={[.16,11]}/><meshStandardMaterial color="#8bb4bf" transparent opacity={.25}/></mesh>
    <mesh position={[6,.2,0]} rotation={[-Math.PI/2,0,Math.PI/2]}><planeGeometry args={[.16,11]}/><meshStandardMaterial color="#d79a83" transparent opacity={.22}/></mesh>
    <mesh position={[0,.47,-1.6]}><boxGeometry args={[1.35,.18,.65]}/><meshStandardMaterial color="#9a8666"/></mesh>
    {trees.map((t,i)=><Tree key={i} {...t} color={t.x<0?"#63a68e":"#ad775f"}/>) }
    <WindTurbine x={-10} z={-3}/><WindTurbine x={-9} z={3.8}/>
    {[0,1,2].map(i=><group key={i} position={[8+i*1.1,.6,-4.4+i*.5]}><mesh><boxGeometry args={[.75,1.2,.75]}/><meshStandardMaterial color={i%2?"#c36b55":"#9a5449"}/></mesh><mesh position={[0,.75,0]}><coneGeometry args={[.55,.5,4]}/><meshStandardMaterial color="#e4b078"/></mesh></group>)}
    <DynamicFrontline x={frontline}/>
  </group>
}
function DynamicFrontline({x}:{x:number}){const ref=useRef<THREE.Group>(null);useFrame(({clock})=>{if(ref.current)ref.current.position.x=THREE.MathUtils.lerp(ref.current.position.x,x,.035);if(ref.current)ref.current.position.y=.32+Math.sin(clock.elapsedTime*2)*.03});return <group ref={ref}><mesh><boxGeometry args={[.12,.08,12]}/><meshStandardMaterial color={GOLD} emissive={GOLD} emissiveIntensity={2}/></mesh><Sparkles count={18} scale={[.6,.7,11]} size={2} color={GOLD} speed={.25}/></group>}

function UnitBody({type,color}:{type:CompanyMarketData["unitType"];color:string}){
  if(type==="aircraft")return <group><mesh rotation={[0,0,Math.PI/2]}><coneGeometry args={[.35,1.45,5]}/><meshStandardMaterial color={color} metalness={.3}/></mesh><mesh><boxGeometry args={[1.4,.12,.28]}/><meshStandardMaterial color="#dce7ee"/></mesh></group>;
  if(type==="drone")return <group>{[-.45,.45].map(x=>[-.35,.35].map(z=><mesh key={`${x}${z}`} position={[x,0,z]}><cylinderGeometry args={[.3,.3,.06,12]}/><meshStandardMaterial color={color}/></mesh>))}<RoundedBox args={[.55,.25,.55]} radius={.12}><meshStandardMaterial color="#dce7ee"/></RoundedBox></group>;
  if(type==="ship")return <group><mesh><boxGeometry args={[1.6,.35,.7]}/><meshStandardMaterial color={color}/></mesh><mesh position={[0,.35,0]}><boxGeometry args={[.55,.38,.45]}/><meshStandardMaterial color="#e1ebee"/></mesh></group>;
  if(type==="command_center")return <group><RoundedBox args={[1.25,.65,1]} radius={.16}><meshStandardMaterial color={color}/></RoundedBox><mesh position={[0,.75,0]}><cylinderGeometry args={[.06,.06,.9,8]}/><meshStandardMaterial color="#d9e9ed"/></mesh><mesh position={[0,.95,0]} rotation={[0,0,.35]}><sphereGeometry args={[.3,10,8,0,Math.PI]}/><meshStandardMaterial color="#a6d4e1"/></mesh></group>;
  return <group><RoundedBox args={[1.35,.42,.85]} radius={.14}><meshStandardMaterial color={color}/></RoundedBox>{type==="tank"&&<><mesh position={[0,.4,0]}><cylinderGeometry args={[.32,.38,.24,10]}/><meshStandardMaterial color={color}/></mesh><mesh position={[-.55,.46,0]} rotation={[0,0,Math.PI/2]}><cylinderGeometry args={[.07,.07,1.4,8]}/><meshStandardMaterial color="#d9e4e8"/></mesh></>}</group>;
}
function CompanyUnit({company,position}:{company:CompanyMarketData;position:[number,number,number]}){
  const {selectedId,range,labels,reducedMotion,set}=useMarketWarStore(); const selected=selectedId===company.id; const group=useRef<THREE.Group>(null); const perf=performance(company,range); const scale=company.status==="private"?.72:unitScale(company.marketCapUsd); const dir=company.country==="USA"?1:-1; const targetX=position[0]+dir*Math.max(-1.4,Math.min(1.4,perf*.18));
  useFrame(({clock})=>{if(!group.current)return; group.current.position.x=THREE.MathUtils.lerp(group.current.position.x,targetX,reducedMotion?.2:.04);group.current.position.y=position[1]+(reducedMotion?0:Math.sin(clock.elapsedTime*1.5+position[2])*.06);group.current.rotation.y=dir===1?Math.PI/2:-Math.PI/2;});
  return <group ref={group} position={position} scale={scale} onClick={e=>{e.stopPropagation();set("selectedId",company.id)}}>
    <mesh position={[0,-.34,0]} rotation={[-Math.PI/2,0,0]}><ringGeometry args={[.64,.8,24]}/><meshBasicMaterial color={selected?GOLD:perf>=0?"#80e7c0":"#ff8a7f"} transparent opacity={selected?1:.45}/></mesh>
    <UnitBody type={company.unitType} color={company.country==="USA"?BLUE:RED}/>
    {perf>2&&<pointLight color="#7cf2c3" intensity={2} distance={2.5}/>} {perf<-1.5&&<Sparkles count={8} scale={1.3} size={2} color="#8b9399" speed={.18}/>} 
    {company.status==="private"&&<mesh position={[0,.05,0]}><boxGeometry args={[1.45,.48,.06]}/><meshBasicMaterial color="#fff" wireframe transparent opacity={.5}/></mesh>}
    {labels&&<Html position={[0,1,0]} center distanceFactor={10} style={{pointerEvents:"none"}}><div className={`unit-label ${selected?"selected":""}`}><span>{company.ticker??"PRIVATE"}</span><b className={perf>=0?"gain":"loss"}>{perf>=0?"+":""}{perf.toFixed(1)}%</b></div></Html>}
  </group>
}
function CameraRig({positions}:{positions:Map<string,[number,number,number]>}){const controls=useRef<any>(null);const selectedId=useMarketWarStore(s=>s.selectedId);const autoCamera=useMarketWarStore(s=>s.autoCamera);const target=useMemo(()=>new THREE.Vector3(),[]);useFrame(()=>{if(!controls.current)return;const p=selectedId?positions.get(selectedId):undefined;target.set(p?.[0]??0,.35,p?.[2]??0);controls.current.target.lerp(target,.045);controls.current.update()});return <OrbitControls ref={controls} makeDefault target={[0,0,0]} minDistance={11} maxDistance={31} minPolarAngle={.4} maxPolarAngle={Math.PI/2.18} enablePan={false} dampingFactor={.06} autoRotate={autoCamera} autoRotateSpeed={.35}/>}
function Scene({companies,usaPercent}:{companies:CompanyMarketData[];usaPercent:number}){const positions=useMemo(()=>unitPositions(companies),[companies]);const shadows=useMarketWarStore(s=>s.shadows);return <><color attach="background" args={["#071118"]}/><fog attach="fog" args={["#071118",21,39]}/><ambientLight intensity={1.7}/><directionalLight position={[-8,15,8]} intensity={2.4} castShadow={shadows}/><BattlefieldTerrain frontline={territoryX(usaPercent)}/>{companies.map(c=><CompanyUnit key={c.id} company={c} position={positions.get(c.id)!}/>)}{shadows&&<ContactShadows position={[0,-.01,0]} opacity={.45} scale={28} blur={2.2}/>}<Environment preset="sunset" environmentIntensity={.35}/><CameraRig positions={positions}/></>}
export function MarketWarScene({companies,usaPercent}:{companies:CompanyMarketData[];usaPercent:number}){return <div className="scene-wrap" aria-hidden="true"><Canvas shadows dpr={[1,1.6]} camera={{position:[15,15,18],fov:42}} gl={{antialias:true,powerPreference:"high-performance"}}><Suspense fallback={null}><Scene companies={companies} usaPercent={usaPercent}/></Suspense></Canvas><div className="scene-vignette"/></div>}
