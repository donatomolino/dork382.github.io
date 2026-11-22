---
title: "The future of tech, from 2027 to 2050"
description: "My vision of the next twenty years: personal AI, temporary software, open protocols, generated interfaces, and smart glasses. What the data supports—and where my own bets begin."
lang: en
slug: the-future-of-tech-2027-2050
translationKey: future-tech
pubDate: 2026-11-22
draft: false
tags:
  - Future
  - Artificial Intelligence
  - Software
  - Interoperability
---

I work with software, digital products, and AI every day. The more time passes, the more I think the biggest change over the next twenty years won't simply be better applications. It will be **our decreasing dependence on applications as we know them today**.

This is my vision, starting in 2027 and extending to more distant scenarios around 2050. I've put data and observable signals alongside each prediction, without treating today's numbers as proof of what will happen tomorrow.

The dates are approximate. Several changes will happen in parallel, and some may never happen at all.

## 2027–2030 Personal AI, not just chatbots

I imagine each of us having one or more personal AIs capable of managing information, devices, services, and everyday tasks. Not just answering questions, but agents that can work with real tools: read a document, compare options, fill out a form, query a service, prepare something, and ask for approval when necessary.

The first signs are already here. In the [McKinsey *State of AI 2025* survey](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), **23%** of respondents said their organizations were scaling at least one agentic AI system; another **39%** reported experimenting with agents. This survey covers organizations, not individuals. It does not mean 62% of people already have personal AI agents.

The ability to operate a computer has also improved dramatically. The [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance) reports a **66.3%** success rate on the OSWorld benchmark in 2025, up from roughly **5%** in 2024. The human reference score is **72.35%**. That's extraordinary progress, but it also shows that an agent may still fail on about one in three tasks in this particular test.

![Chart comparing AI performance on the OSWorld benchmark: roughly 5% in 2024, 66.3% in 2025, and a human reference of 72.35%.](/blog/charts/en/agents-osworld.svg)

*Chart 1 — Data: [Stanford HAI, AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance). This benchmark measures specific computer tasks, not the overall reliability of every AI agent.*

For me, the most interesting part is independence from providers. I want an AI that keeps its memory and context even when I switch models, and that can operate locally or in the cloud. Open standards are moving in a promising direction: in December 2025, the [Model Context Protocol joined the Agentic AI Foundation](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), under the Linux Foundation.

But none of this is solved. Portable memory, permissions, security, accountability, and cost remain open problems. In 2025, [Gartner predicted](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) that more than 40% of agentic AI projects would be canceled by the end of 2027. **That's Gartner's forecast, not a failure rate that has already been observed.**

## 2028–2033 More personalized, sometimes temporary software

I don't think Salesforce, SAP, ServiceNow, or Microsoft Dynamics will disappear. I think their role will change.

Today, we often adapt how we work to the screens of the software we bought. I imagine the opposite: software adapting to the workflow, built or modified quickly when necessary. An internal application might last a few months, be replaced by something else, and leave the data, APIs, identities, and core business logic intact.

Code generation is reducing the cost of building interfaces and automations. In 2025, Gartner [forecast](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) that 33% of enterprise software applications would include agentic AI by 2028, up from less than 1% in 2024. That's still a **projection**, not a measurement of how long applications actually last.

My bet isn't that everything will be constantly rewritten. That would be a maintenance nightmare. It's that **the visible layer of software will become cheaper to change**, while data contracts, permissions, integrations, controls, and infrastructure will need to become much more stable. In regulated industries, audits, security, and reliability will still limit how quickly things can change.

## 2028–2033 Phones, tablets, and PCs converge

I think that for most everyday tasks, the smartphone will eventually become the only computer we need: a phone in our pocket, a tablet when unfolded, and a full desktop when connected to a monitor.

This isn't science fiction. In March 2026, [Google documented](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays) Android desktop experiences on external monitors, developed in collaboration with Samsung. The problem is no longer demonstrating that a phone can show desktop windows. It's making the whole experience reliable enough to replace a laptop.

Real obstacles remain: professional applications, peripherals, battery life, heat dissipation, heavy multitasking, hardware architectures, and compatibility. Game developers, 3D artists, and people working in specialized environments may still need dedicated workstations.

I don't believe the PC will die. **I believe that, for many people, the PC will become a mode of use rather than a separate device.**

## 2030–2037 Shared protocols instead of closed ecosystems

This is probably the prediction I care about most.

I should be able to use one messaging app while somebody else uses a completely different client, and we should still be able to communicate. Just like we already do with email. I can imagine the same principle gradually extending to social networks, content, digital services, and AI.

Some pieces already exist. [ActivityPub](https://www.w3.org/TR/activitypub/) has been a W3C recommendation for decentralized social networking since 2018. The EU's [Digital Markets Act](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en) requires designated large messaging platforms to provide specific forms of interoperability; the first third-party services interoperable with WhatsApp arrived in November 2025. In the agent ecosystem, standards like MCP address how tools and data can connect.

These are signals, not guarantees. Interoperability raises challenges around moderation, spam, encryption, identity, and business incentives. Platforms have powerful commercial reasons to keep people inside their ecosystems.

Still, I believe value will gradually move **from having to use the same app toward being free to choose the interface you prefer**.

## 2030–2038 Interfaces generated for the task at hand

Today, we use products made of predefined screens that are more or less identical for everyone. Tomorrow, we could describe the result we want and receive an interface built specifically around that task: just the relevant data, the necessary actions, and controls appropriate for the device.

We are not starting from zero. [Google Research has explored dynamically generated interfaces](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), while the open [A2UI project](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/) defines a way to describe agent-produced interfaces and render them across different frameworks.

Imagine opening my AI and asking it to compare household expenses. Instead of navigating three different apps, it generates a chart, a table, and two useful actions. When I change my goal, the interface changes too.

This doesn't mean that every pixel should be improvised. Sensitive operations will still need dependable components, accessibility, consistency, and explicit confirmation. **An interface can be dynamic without being unpredictable.**

## 2030–2040 Sensors almost everywhere

Glasses, headphones, rings, clothing, and household objects will collect increasing amounts of information about the environment and our activities. My personal AI could combine some of these signals to understand far more context than a simple written prompt provides.

One number puts the scale into perspective. According to [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/), there were **18.5 billion** connected IoT devices worldwide in 2024; the estimate for 2025 was **21.1 billion**, and the forecast for 2030 reaches **39 billion**.

![Chart showing global IoT connections: 18.5 billion in 2024, an estimated 21.1 billion in 2025, and a forecast 39 billion in 2030.](/blog/charts/en/iot-devices.svg)

*Chart 2 — [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/). These totals include industry, infrastructure, and other connected devices—not just wearables. The 2030 figure is a forecast.*

But sensors do not automatically produce useful insights. Inaccurate or out-of-context data can be worse than no data at all. And when a system knows your habits, location, and personal signals, privacy and security cannot be afterthoughts.

The future I want involves **granular control, local processing whenever possible, and minimal data sharing**. Not permanent surveillance disguised as convenience.

## 2032–2040 Digital identity independent of platforms

Today, our accounts, preferences, files, and histories are scattered across services that often don't communicate. I imagine a clearer separation between **who I am**, **which data I control**, and **which application I use to access it**.

Europe is already moving in this direction. The [European Digital Identity Wallet](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en) is intended to be made available by member states by the end of 2026, offering ways to identify ourselves and share specific credentials. The [Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act), applicable since September 2025, also introduces measures for accessing data from connected products and switching cloud providers.

There is an important distinction: **an identity wallet is not the same thing as a universal profile containing every conversation, preference, and AI memory**. That broader kind of portability remains to be built, both technically and legally.

My prediction is that the center of gravity will gradually shift from platform accounts toward user control. It won't take just one standard, and it won't happen everywhere at the same speed.

## 2035–2045 Truly spatial augmented reality

I don't think the endgame for AR is having 2D notifications and windows floating in front of our eyes. I expect digital objects that understand the geometry of the world, stay anchored in place, and respond to physical spaces.

Think of instructions attached to something I'm repairing, information about an object I'm looking at, or virtual work surfaces placed wherever I need them. That's very different from projecting a tiny screen into your field of vision.

Device sales shouldn't be confused with the maturity of that vision. [IDC](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/) forecasts around **3 million units shipped in 2026** for glasses with displays, rising to **12.2 million in 2030**. That category also includes devices that put images in front of your eyes without delivering full spatial AR.

The real bottlenecks will still be battery life, heat, weight, field of view, precise tracking, optical quality, and social acceptance. That's why I put truly mature AR further into the future.

## 2040–2050+ Glasses could replace the smartphone

This is my most distant prediction, and the one where I'm least confident about the timing.

At some point, I think glasses could take over many of the phone's functions: communication, navigation, contextual information, translation, photography, and interaction with a personal AI. We might no longer need to pull out a screen for every small action.

The market is moving, although we're still far from a replacement. IDC forecasts that **display-less smart glasses** will grow from **13.6 million annual shipments in 2026** to **27.3 million in 2030**. Glasses with displays are a separate, smaller market.

![Chart showing IDC's forecast smart-glasses shipments with and without displays in 2026 and 2030.](/blog/charts/en/smart-glasses.svg)

*Chart 3 — [IDC, 2026](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/). These are forecasts for annual shipments, not active users or the percentage of smartphones replaced.*

The numbers show market interest and projected growth. **They do not prove that smartphones will disappear.** For extended reading, complex work, or a private display, we may still prefer other devices.

My bet is more specific: glasses could become the primary interface for many brief, context-aware interactions. The smartphone could remain, but move into the background.

## The real shift: from apps to connected systems

When I put these nine predictions together, I see a common thread:

- **Personal AIs** could become the starting point for many tasks.
- **Data, identity, and protocols** should remain usable when we change software.
- **Interfaces** could adapt to the job instead of forcing us through fixed screens.
- **Devices** could become different ways of accessing the same digital environment.

I don't believe in a future without apps. I believe in one where apps are **less central and easier to replace**.

Of course, the opposite could happen: closed personal AIs, even more concentrated data, invasive sensors, and new dependencies on a handful of providers. Technology alone doesn't guarantee openness. That takes standards, regulation, well-designed products, and people demanding control.

That's why the most important prediction here isn't about a device or a date.

**The real value won't lie only in the software we use, but in the ability of our systems to communicate without taking away our freedom to choose.**

---

### Sources and methodology

I used information published by **November 22, 2026**. Each chart cites its underlying source. Estimates and market forecasts are labeled as such; the date ranges in this essay are **my personal hypotheses**, not statistical projections.

Main sources: [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report), [McKinsey — State of AI 2025](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), [Gartner — Agentic AI](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027), [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/), [IDC — Smart Glasses](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/), [Android Developers](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays), [W3C ActivityPub](https://www.w3.org/TR/activitypub/), [Linux Foundation — Agentic AI Foundation](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), [Google Research — Generative UI](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), [European Commission — Digital Identity](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en), [European Commission — Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act), and [Digital Markets Act — Interoperability](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en).
