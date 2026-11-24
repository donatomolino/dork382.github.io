---
title: "The future of tech, from 2027 to 2050"
description: "My vision of the next twenty years: personal AI, temporary software, open protocols, generated interfaces, and smart glasses. What the data supports—and where my own bets begin."
lang: en
slug: the-future-of-tech-2027-2050
translationKey: future-tech
pubDate: 2026-11-22
updatedDate: 2026-11-24
draft: false
sources:
  - https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance
  - https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025
  - https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027
  - https://www.gartner.com/en/newsroom/press-releases/2026-07-01-gartner-says-us-dollars-234-billion-in-enterprise-application-software-spend-is-at-risk-from-agentic-artificial-intelligence
  - https://iot-analytics.com/number-connected-iot-devices/
  - https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/
  - https://www.w3.org/TR/activitypub/
  - https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation
  - https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en
  - https://digital-strategy.ec.europa.eu/en/policies/data-act
  - https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en
  - https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays
  - https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/
  - https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/
  - https://help.openai.com/en/articles/11752874-chatgpt-agent
  - https://help.openai.com/en/articles/11487775-connectors-in-chatgpt
  - https://www.anthropic.com/research/claude-code-best-practices
  - https://www.anthropic.com/news/model-context-protocol
  - https://help.openai.com/en/articles/12584461-developer-mode-and-full-mcp-connectors-in-chatgpt-beta
  - https://docs.github.com/en/copilot/get-started/about-github-copilot
  - https://www.samsung.com/it/apps/samsung-dex/
  - https://www.home-assistant.io/integrations/matter
  - https://developer.apple.com/visionos/
  - https://www.meta.com/design-at-meta/blog/three-design-leaders-share-their-vision-for-wearables-meta/
  - https://www.meta.com/ai-glasses/
  - https://www.meta.com/ai-glasses/meta-ray-ban-display/
  - https://developers.meta.com/blog/meta-connect-recap-ai-glasses/
  - https://www.gartner.com/en/newsroom/press-releases/2025-08-26-gartner-predicts-40-percent-of-enterprise-apps-will-feature-task-specific-ai-agents-by-2026-up-from-less-than-5-percent-in-2025
tags:
  - Future
  - Artificial Intelligence
  - Software
  - Interoperability
---

I work with software, digital products, and AI every day. The more time passes, the more I think the biggest change over the next twenty years won't simply be better applications. It will be **our decreasing dependence on applications as we know them today**.

This is my vision, starting in 2027 and extending to more distant scenarios around 2050. I've put data and observable signals alongside each prediction, without treating today's numbers as proof of what will happen tomorrow.

The dates are approximate. Several changes will happen in parallel, and some may never happen at all.

**How to read this essay:** the charts distinguish measured results, survey responses, analyst forecasts, and my own conceptual diagrams. A forecast is not a fact about the future, and a diagram is not a statistic. The examples below are either documented current capabilities or explicitly hypothetical product scenarios.

## 2027–2030 Personal AI, not just chatbots

I imagine each of us having one or more personal AIs capable of managing information, devices, services, and everyday tasks. Not just answering questions, but agents that can work with real tools: read a document, compare options, fill out a form, query a service, prepare something, and ask for approval when necessary.

The products I watch today are [ChatGPT](https://help.openai.com/en/articles/11752874-chatgpt-agent) and Claude. ChatGPT can use [connected apps](https://help.openai.com/en/articles/11487775-connectors-in-chatgpt) to retrieve information from services such as Google Drive or Slack; Claude Code can work on a repository and use [MCP tools](https://www.anthropic.com/research/claude-code-best-practices) with permissions. Those are concrete examples of the move from conversation to action. They still differ in access, memory and capabilities, and they are not interchangeable personal operating systems.

The first signs are already here. In the [McKinsey *State of AI 2025* survey](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), **23%** of respondents said their organizations were scaling at least one agentic AI system; another **39%** reported experimenting with agents. This survey covers organizations, not individuals. It does not mean 62% of people already have personal AI agents.

![Two separate 2025 survey figures: 23% of respondents report their organizations scaling agentic AI, and an additional 39% report experimenting.](/blog/charts/en/enterprise-agents-survey.svg)

*Chart 1 — Organizational adoption signals, not consumer adoption. Source: [McKinsey, State of AI 2025](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), global survey of 1,993 respondents, June–July 2025. The groups describe reported stages; they are not a forecast.*

The ability to operate a computer has also improved dramatically. The [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance) reports a **66.3%** success rate on the OSWorld benchmark in 2025, up from roughly **5%** in 2024. The human reference score is **72.35%**. That's extraordinary progress, but it also shows that an agent may still fail on about one in three tasks in this particular test.

![Chart comparing AI performance on the OSWorld benchmark: roughly 5% in 2024, 66.3% in 2025, and a human reference of 72.35%.](/blog/charts/en/agents-osworld.svg)

*Chart 2 — Data: [Stanford HAI, AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance). This benchmark measures specific computer tasks, not the overall reliability of every AI agent.*

**A concrete example:** I ask my personal AI to plan a weekend trip. It reads a train timetable, compares hotels against my budget, proposes two itineraries, and prepares a reservation. Crucially, it **doesn't pay or book until I approve**. This is a hypothetical end-to-end scenario, not a claim that today's agents can do it reliably.

For me, the most interesting part is independence from providers. I want an AI that keeps its memory and context even when I switch models, and that can operate locally or in the cloud. Open standards are moving in a promising direction: in December 2025, the [Model Context Protocol joined the Agentic AI Foundation](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), under the Linux Foundation.

But none of this is solved. Portable memory, permissions, security, accountability, and cost remain open problems. In 2025, [Gartner predicted](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) that more than 40% of agentic AI projects would be canceled by the end of 2027. **That's Gartner's forecast, not a failure rate that has already been observed.**

## 2028–2033 More personalized, sometimes temporary software

I don't think Salesforce, SAP, ServiceNow, or Microsoft Dynamics will disappear. I think their role will change.

[Claude Code](https://www.anthropic.com/research/claude-code-best-practices) and [GitHub Copilot](https://docs.github.com/en/copilot/get-started/about-github-copilot) already let developers delegate parts of a coding workflow to agents. Copilot, for example, can explore a repository, make changes and prepare a pull request for review. I see these products as the beginning of cheaper, more adaptable software, not as proof that companies can dispense with engineering, testing or maintenance.

Today, we often adapt how we work to the screens of the software we bought. I imagine the opposite: software adapting to the workflow, built or modified quickly when necessary. An internal application might last a few months, be replaced by something else, and leave the data, APIs, identities, and core business logic intact.

Code generation is reducing the cost of building interfaces and automations. In 2025, Gartner [forecast](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) that 33% of enterprise software applications would include agentic AI by 2028, up from less than 1% in 2024. That's still a **projection**, not a measurement of how long applications actually last.

![Two separate Gartner indicators, 2024–2028: agentic AI, below 1% in 2024 and 33% forecast for 2028; task-specific agents, below 5% in 2025 and up to 40% forecast by end-2026. Missing comparable annual figures marked N/A.](/blog/charts/en/enterprise-agent-forecast.svg)

*Chart 3 — Two distinct measures: [Gartner, June 2025](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) for agentic AI in applications and [Gartner, August 2025](https://www.gartner.com/en/newsroom/press-releases/2025-08-26-gartner-predicts-40-percent-of-enterprise-apps-will-feature-task-specific-ai-agents-by-2026-up-from-less-than-5-percent-in-2025) for integrated task-specific agents. The 2026 and 2028 figures are forecasts issued in 2025, not measured outcomes. N/A indicates that I found no comparable annual figure, not zero.*

For 2027, Gartner also predicts that roughly one-third of agentic AI implementations will combine agents with different skills. That is a third measure, concerning implementations rather than the share of applications, so I do not plot it as part of either series.

My bet isn't that everything will be constantly rewritten. That would be a maintenance nightmare. It's that **the visible layer of software will become cheaper to change**, while data contracts, permissions, integrations, controls, and infrastructure will need to become much more stable. In regulated industries, audits, security, and reliability will still limit how quickly things can change.

Here's how I picture it in commerce: a retailer creates a temporary dashboard for a two-week campaign. The interface is tailored to the team running the launch, but inventory, customer records, payment permissions, and order APIs remain in established systems. When the campaign ends, the dashboard can go away **without taking the data with it**.

This isn't just about developer productivity. In July 2026, [Gartner estimated](https://www.gartner.com/en/newsroom/press-releases/2026-07-01-gartner-says-us-dollars-234-billion-in-enterprise-application-software-spend-is-at-risk-from-agentic-artificial-intelligence) that **up to $234 billion in enterprise application SaaS spending**, approximately **20%** of that category, could be exposed to agentic arbitrage by 2030: agents doing work across tools without users opening each interface. That money is **exposed to a different way of buying software**, not already lost revenue.

![Diagram of Gartner's 2030 enterprise application SaaS exposure estimate: up to $234 billion, roughly 20% of the category.](/blog/charts/en/saas-exposure-2030.svg)

*Chart 4 — Source: [Gartner, July 2026](https://www.gartner.com/en/newsroom/press-releases/2026-07-01-gartner-says-us-dollars-234-billion-in-enterprise-application-software-spend-is-at-risk-from-agentic-artificial-intelligence). Market exposure is an analyst estimate about business-model pressure, **not a prediction of $234B in losses**.*

### A new market: tools for AI

I believe a market will emerge in which the customer is still the person, but their AI increasingly operates the software on their behalf. If ChatGPT, Claude or a personal assistant can select and combine the tools it needs, building an entire application merely to make the user navigate its interface becomes less compelling.

Think of a document-conversion service, a specialized search engine, a real-time price-comparison service, or a tool that generates invoices and checks payments. An agent could invoke them through APIs or MCP and assemble the experience it needs. I wouldn't necessarily buy four apps: I'd pay, directly or through my assistant, for the capabilities I actually need.

**The product would no longer necessarily be the application, but what it lets you do.** Value would lie in the data, output quality, reliability, permissions and interoperability. I can imagine usage-based, per-operation or outcome-based business models alongside traditional subscriptions. That doesn't mean complete products or interfaces will become pointless; I am thinking particularly of specialized services that an agent can compose.

It's an idea I want to explore in more depth in a future article.

## 2028–2033 Phones, tablets, and PCs converge

I think that for most everyday tasks, the smartphone will eventually become the only computer we need: a phone in our pocket, a tablet when unfolded, and a full desktop when connected to a monitor.

[Samsung DeX](https://www.samsung.com/it/apps/samsung-dex/) is one of the clearest examples: compatible Galaxy phones, including foldable models, can offer a desktop-style workspace on an external screen. A foldable phone with that capability is close to my idea of a pocket computer. The remaining question is whether the software ecosystem can support all the work people currently do on a PC.

This isn't science fiction. In March 2026, [Google documented](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays) Android desktop experiences on external monitors, developed in collaboration with Samsung. The problem is no longer demonstrating that a phone can show desktop windows. It's making the whole experience reliable enough to replace a laptop.

**What I mean in practice:** I start editing a document on my phone during a train ride. At home I connect it to a monitor and keyboard, and the same device opens a multi-window workspace. Google and Samsung's [Android 16 connected-display work](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays) is a real technical step in this direction. That doesn't mean every Android phone supports the feature or that every desktop application runs there.

Real obstacles remain: professional applications, peripherals, battery life, heat dissipation, heavy multitasking, hardware architectures, and compatibility. Game developers, 3D artists, and people working in specialized environments may still need dedicated workstations.

I don't believe the PC will die. **I believe that, for many people, the PC will become a mode of use rather than a separate device.**

## 2030–2037 Shared protocols instead of closed ecosystems

This is probably the prediction I care about most.

The [Model Context Protocol (MCP)](https://www.anthropic.com/news/model-context-protocol), originally introduced by Anthropic, is particularly interesting to me. Imagine a GitHub tool or a private document service exposed through an MCP server: compatible assistants can discover and invoke the operations that server exposes, within its permissions. [ChatGPT supports MCP-based apps](https://help.openai.com/en/articles/12584461-developer-mode-and-full-mcp-connectors-in-chatgpt-beta), and Claude Code can connect to MCP servers. The protocol does not make their models, memory or identities portable; it standardizes a connection to tools and context. I would like that same separation to spread much further.

I should be able to use one messaging app while somebody else uses a completely different client, and we should still be able to communicate. Just like we already do with email. I can imagine the same principle gradually extending to social networks, content, digital services, and AI.

Some pieces already exist. [ActivityPub](https://www.w3.org/TR/activitypub/) has been a W3C recommendation for decentralized social networking since 2018. The EU's [Digital Markets Act](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en) requires designated large messaging platforms to provide specific forms of interoperability; the first third-party services interoperable with WhatsApp arrived in November 2025. In the agent ecosystem, standards like MCP address how tools and data can connect.

![Timeline of interoperability milestones: ActivityPub in 2018, EU Data Act and third-party WhatsApp chats in 2025, MCP foundation in 2025, and EU identity wallets due by end of 2026.](/blog/charts/en/interoperability-milestones.svg)

*Chart 5 — Selected, independently sourced milestones: [W3C ActivityPub](https://www.w3.org/TR/activitypub/), [EU Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act), [EU messaging interoperability](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en), [Linux Foundation / MCP](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), and [EU Digital Identity](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en). The last item is a regulatory availability deadline, not evidence of universal adoption.*

**A concrete distinction:** two people can follow each other through compatible ActivityPub-powered social services, such as Mastodon servers, even when they use different services. That's not the same as being able to message freely between every popular chat app. The [EU's November 2025 WhatsApp integrations](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en) are a narrower, regulated step, with important limitations.

These are signals, not guarantees. Interoperability raises challenges around moderation, spam, encryption, identity, and business incentives. Platforms have powerful commercial reasons to keep people inside their ecosystems.

Still, I believe value will gradually move **from having to use the same app toward being free to choose the interface you prefer**.

## 2030–2038 Interfaces generated for the task at hand

Today, we use products made of predefined screens that are more or less identical for everyone. Tomorrow, we could describe the result we want and receive an interface built specifically around that task: just the relevant data, the necessary actions, and controls appropriate for the device.

[Apps in ChatGPT](https://help.openai.com/en/articles/11487775-connectors-in-chatgpt) are an early consumer-facing example of bringing third-party services into a conversation. The next step I imagine is an interface assembled from safe components for each task, rather than merely embedding an existing app. Google's A2UI is particularly relevant to that distinction.

We are not starting from zero. [Google Research has explored dynamically generated interfaces](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), while the open [A2UI project](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/) defines a way to describe agent-produced interfaces and render them across different frameworks.

![Concept diagram showing a generated-interface workflow: user request, permissioned data, approved components, and explicit human confirmation.](/blog/charts/en/generated-interface-flow.svg)

*Chart 6 — **Conceptual design**, not measured adoption or a live product. Informed by [Google Research's generative UI work](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/) and the [A2UI project](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/).*

Imagine opening my AI and asking it to compare household expenses. Instead of navigating three different apps, it generates a chart, a table, and two useful actions. More specifically, I could ask: *“Why was my electricity bill higher this month?”* It would request permission to read the relevant bills and meter data, show a breakdown of consumption versus tariffs, and offer a comparison with last month. Any payment or supplier change would still require an explicit decision from me.

Google's own [A2UI examples](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/) show agents returning structured interactive components instead of plain chat text. That's a real technical direction, although secure cross-provider execution remains harder than displaying a useful chart. When I change my goal, the interface changes too.

This doesn't mean that every pixel should be improvised. Sensitive operations will still need dependable components, accessibility, consistency, and explicit confirmation. **An interface can be dynamic without being unpredictable.**

## 2030–2040 Sensors almost everywhere

Glasses, headphones, rings, clothing, and household objects will collect increasing amounts of information about the environment and our activities. My personal AI could combine some of these signals to understand far more context than a simple written prompt provides.

A concrete example is [Home Assistant with Matter](https://www.home-assistant.io/integrations/matter): compatible devices can communicate locally, without needing every manufacturer’s separate cloud integration. Combine that kind of home infrastructure with wearable devices such as [Ray-Ban Meta](https://www.meta.com/ai-glasses/), and I can imagine an assistant with richer situational context. Combining these systems securely is my proposed future, not a feature they already provide together.

One number puts the scale into perspective. According to [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/), there were **18.5 billion** connected IoT devices worldwide in 2024; the estimate for 2025 was **21.1 billion**, and the forecast for 2030 reaches **39 billion**.

![Chart showing global IoT connections: 18.5 billion in 2024, an estimated 21.1 billion in 2025, and a forecast 39 billion in 2030.](/blog/charts/en/iot-devices.svg)

*Chart 7 — [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/). These totals include industry, infrastructure, and other connected devices—not just wearables. The 2030 figure is a forecast.*

A useful home example doesn't need to be futuristic: a door sensor, a temperature sensor, and a small local server can already provide a simple activity history and automations. My longer-term bet is that a personal AI could connect those signals to genuinely useful tasks **without sending an entire household's raw data to a third party**. That's a design goal, not something to assume by default.

But sensors do not automatically produce useful insights. Inaccurate or out-of-context data can be worse than no data at all. And when a system knows your habits, location, and personal signals, privacy and security cannot be afterthoughts.

The future I want involves **granular control, local processing whenever possible, and minimal data sharing**. Not permanent surveillance disguised as convenience.

## 2032–2040 Digital identity independent of platforms

Today, our accounts, preferences, files, and histories are scattered across services that often don't communicate. I imagine a clearer separation between **who I am**, **which data I control**, and **which application I use to access it**.

Europe is already moving in this direction. The [European Digital Identity Wallet](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en) is intended to be made available by member states by the end of 2026, offering ways to identify ourselves and share specific credentials. The [Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act), applicable since September 2025, also introduces measures for accessing data from connected products and switching cloud providers.

Consider a practical use case: proving a specific qualification or eligibility to a service without creating another full account. This is the narrower problem digital credentials can help with. It is very different from handing over a complete profile of who I am.

There is an important distinction: **an identity wallet is not the same thing as a universal profile containing every conversation, preference, and AI memory**. That broader kind of portability remains to be built, both technically and legally.

My prediction is that the center of gravity will gradually shift from platform accounts toward user control. It won't take just one standard, and it won't happen everywhere at the same speed.

## 2035–2045 Truly spatial augmented reality

I don't think the endgame for AR is having 2D notifications and windows floating in front of our eyes. I expect digital objects that understand the geometry of the world, stay anchored in place, and respond to physical spaces.

[Apple Vision Pro and visionOS](https://developer.apple.com/visionos/) already demonstrate spatial applications, while [Meta Orion](https://www.meta.com/design-at-meta/blog/three-design-leaders-share-their-vision-for-wearables-meta/) was presented as a prototype of AR glasses. These are two distinct approaches: a headset available as a product and a prototype exploring a glasses form factor. Neither demonstrates that lightweight, all-day spatial AR is a solved engineering problem.

Picture fixing a machine: instead of looking down at a video, I see the correct component highlighted in its actual physical location, a torque value beside the fastener, and the next step only after confirming the previous one. **That's an illustrative future workflow**, not something I assume consumer AR glasses can deliver reliably today.

Think of instructions attached to something I'm repairing, information about an object I'm looking at, or virtual work surfaces placed wherever I need them. That's very different from projecting a tiny screen into your field of vision.

Device sales shouldn't be confused with the maturity of that vision. [IDC](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/) forecasts around **3 million units shipped in 2026** for glasses with displays, rising to **12.2 million in 2030**. That category also includes devices that put images in front of your eyes without delivering full spatial AR.

The real bottlenecks will still be battery life, heat, weight, field of view, precise tracking, optical quality, and social acceptance. That's why I put truly mature AR further into the future.

## 2040–2050+ Glasses could replace the smartphone

This is my most distant prediction, and the one where I'm least confident about the timing.

The difference between [Ray-Ban Meta](https://www.meta.com/ai-glasses/) and [Meta Ray-Ban Display](https://www.meta.com/ai-glasses/meta-ray-ban-display/) matters here. The former illustrates voice-first AI, cameras and hands-free capture without a display; the latter adds an in-lens display and the Meta Neural Band. Meta has also [documented web-app experiences for its display glasses](https://developers.meta.com/blog/meta-connect-recap-ai-glasses/). I see these products as early steps toward a new interface, not a complete replacement for the smartphone. A small display in front of the eye is also different from spatial AR that anchors content throughout a room.

At some point, I think glasses could take over many of the phone's functions: communication, navigation, contextual information, translation, photography, and interaction with a personal AI. We might no longer need to pull out a screen for every small action.

The market is moving, although we're still far from a replacement. IDC forecasts that **display-less smart glasses** will grow from **13.6 million annual shipments in 2026** to **27.3 million in 2030**. Glasses with displays are a separate, smaller market.

![Chart showing IDC's forecast smart-glasses shipments with and without displays in 2026 and 2030.](/blog/charts/en/smart-glasses.svg)

*Chart 8 — [IDC, 2026](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/). These are forecasts for annual shipments, not active users or the percentage of smartphones replaced.*

The use case I find most convincing is small, frequent assistance: walking directions, a translation while traveling, or a reminder linked to a physical place, without taking out a phone. The harder question isn't whether a device can display that information. It's whether I would trust it to be comfortable, private, reliable, and socially acceptable all day.

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

Imagine the full path: I ask a personal AI to organize my next trip. It accesses travel services through authorized connections, presents a temporary interface to compare routes, uses a portable credential only where necessary, and later gives me directions on whichever screen I have available. **The experience is continuous even if the apps, services, and devices behind it are not.** That's the future I'm betting on, not one magic product replacing everything.

That's why the most important prediction here isn't about a device or a date.

**The real value won't lie only in the software we use, but in the ability of our systems to communicate without taking away our freedom to choose.**

---

### Sources and methodology

This revision includes sources available by **November 24, 2026**. Eight charts are embedded as small static SVGs, with a direct source and a scope note beneath each. Charts 1, 2, 3, 4, 7 and 8 show reported data, survey findings or cited analyst estimates. Chart 5 is a timeline of sourced milestones; Chart 6 is a **conceptual product workflow**, not measured data.

The 2027–2050 dates are **my hypotheses**, not statistical forecasts. Enterprise surveys do not measure consumer behavior; shipment forecasts do not equal active users; and dollars **exposed** to disruption are not dollars **lost**. I deliberately keep those distinctions visible.

Primary sources by subject:

- **AI performance and enterprise adoption:** [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance), [McKinsey State of AI 2025](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), [Gartner on agents in software (2025)](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027).
- **Software economics and interfaces:** [Gartner's July 2026 SaaS exposure estimate](https://www.gartner.com/en/newsroom/press-releases/2026-07-01-gartner-says-us-dollars-234-billion-in-enterprise-application-software-spend-is-at-risk-from-agentic-artificial-intelligence), [Android Developers on connected displays](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays), [Google Research on generative UI](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), [Google Developers on A2UI](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/).
- **Protocols, digital identity and policy:** [W3C ActivityPub](https://www.w3.org/TR/activitypub/), [Linux Foundation / MCP](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), [EU Digital Markets Act interoperability](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en), [European Digital Identity Wallet](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en), [EU Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act).
- **Connected devices and glasses:** [IoT Analytics](https://iot-analytics.com/number-connected-iot-devices/), [IDC smart glasses market outlook (2026)](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/).
