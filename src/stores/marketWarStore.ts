import { create } from "zustand";
import { persist } from "zustand/middleware";
import type { Country, TimeRange, UnitType } from "../features/market-war/types";

type Quality="low"|"high";
interface State { selectedId?:string; range:TimeRange; count:10|25|50|100|"all"; search:string; country:"ALL"|Country; sector:string; unit:"ALL"|UnitType; sign:"ALL"|"POSITIVE"|"NEGATIVE"; showPrivate:boolean; settingsOpen:boolean; filtersOpen:boolean; tableOpen:boolean; quality:Quality; particles:boolean; shadows:boolean; labels:boolean; autoCamera:boolean; reducedMotion:boolean; demoMotion:boolean; set:<K extends keyof State>(key:K,value:State[K])=>void; }
export const useMarketWarStore=create<State>()(persist((set)=>({range:"24H",count:25,search:"",country:"ALL",sector:"ALL",unit:"ALL",sign:"ALL",showPrivate:true,settingsOpen:false,filtersOpen:false,tableOpen:false,quality:"high",particles:true,shadows:true,labels:true,autoCamera:false,reducedMotion:false,demoMotion:false,set:(key,value)=>set({[key]:value})}),{name:"market-war-settings",partialize:s=>({quality:s.quality,particles:s.particles,shadows:s.shadows,labels:s.labels,autoCamera:s.autoCamera,reducedMotion:s.reducedMotion,showPrivate:s.showPrivate})}));
