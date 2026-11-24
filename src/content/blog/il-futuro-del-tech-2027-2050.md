---
title: "Il futuro del tech, dal 2027 al 2050"
description: "La mia visione dei prossimi vent'anni: AI personali, software temporaneo, protocolli aperti, interfacce generate e smart glasses. Cosa dicono i dati, e dove iniziano le mie scommesse."
lang: it
slug: il-futuro-del-tech-2027-2050
translationKey: future-tech
pubDate: 2026-11-22
updatedDate: 2026-11-24
draft: false
sources:
  - https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance
  - https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025
  - https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027
  - https://iot-analytics.com/Number-Connected-Iot-Devices/
  - https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/
  - https://www.w3.org/TR/activitypub/
  - https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation
  - https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en
  - https://digital-strategy.ec.europa.eu/en/policies/data-act
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
tags:
  - Futuro
  - Intelligenza artificiale
  - Software
  - Interoperabilità
---

Lavoro ogni giorno con software, prodotti digitali e AI. Più passa il tempo, più penso che il cambiamento dei prossimi vent'anni non sarà semplicemente avere applicazioni migliori: sarà **dipendere meno dalle applicazioni come le conosciamo oggi**.

Questa è la mia visione, dal 2027 fino agli scenari più lontani del 2050. Ho cercato di mettere accanto alle previsioni dati e segnali che possiamo già osservare nel 2026, senza confondere i numeri di oggi con certezze sul futuro.

Le date sono indicative. Molte evoluzioni avverranno in parallelo e alcune potrebbero non avvenire affatto.

## 2027–2030 Avremo AI personali, non soltanto chatbot

Immagino che ognuno di noi avrà una o più AI personali capaci di gestire informazioni, dispositivi, servizi e attività. Non soltanto risposte a domande, ma agenti che lavorano con strumenti reali: leggono un documento, confrontano opzioni, compilano una richiesta, interrogano un servizio, preparano qualcosa e chiedono conferma quando serve.

I prodotti che osservo oggi sono soprattutto [ChatGPT](https://help.openai.com/en/articles/11752874-chatgpt-agent) e Claude. ChatGPT può utilizzare [app collegate](https://help.openai.com/en/articles/11487775-connectors-in-chatgpt) per recuperare informazioni da servizi come Google Drive o Slack; Claude Code può lavorare su un repository e usare [strumenti MCP](https://www.anthropic.com/research/claude-code-best-practices) con permessi specifici. Sono esempi concreti del passaggio dalla conversazione all'azione. Non sono ancora sistemi operativi personali intercambiabili: cambiano memoria, accessi e capacità.

Il segnale è già visibile. Nel sondaggio [McKinsey *State of AI 2025*](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), il **23%** dei rispondenti indicava che la propria organizzazione stava ampliando l'uso di almeno un sistema AI agentico; un ulteriore **39%** dichiarava sperimentazioni. È un'indagine sulle organizzazioni, non sulle persone: non significa che il 62% della popolazione abbia già un agente personale.

Anche la capacità di operare su un computer è cambiata rapidamente. Lo [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance) riporta un successo del **66,3%** sul benchmark OSWorld nel 2025, contro circa il **5%** nel 2024. Il riferimento umano è **72,35%**. È una progressione notevole, ma significa anche che un agente può ancora fallire circa una volta su tre in questo tipo di test.

![Grafico: successo degli agenti AI su OSWorld, circa 5% nel 2024, 66,3% nel 2025 e 72,35% per il riferimento umano.](/blog/charts/agents-osworld.svg)

*Grafico 1 — Dati: [Stanford HAI, AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report/technical-performance). Il benchmark misura compiti specifici al computer; non dimostra l'affidabilità degli agenti in ogni situazione.*

Per me la parte più interessante è l'indipendenza dai fornitori. Vorrei un'AI che possa conservare memoria e contesto anche quando cambio modello, e che possa usare strumenti locali oppure cloud. La direzione degli standard è promettente: nel dicembre 2025 il [Model Context Protocol è entrato nella Agentic AI Foundation](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), sotto la Linux Foundation.

Ma non è tutto risolto. Memoria portabile, permessi, sicurezza, responsabilità e costi rimangono problemi aperti. E Gartner [prevedeva nel 2025](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) che oltre il 40% dei progetti agentici sarebbe stato abbandonato entro il 2027. **È una previsione di Gartner, non un tasso di fallimento già osservato.**

## 2028–2033 Software più personalizzato, e a volte temporaneo

Non penso che Salesforce, SAP, ServiceNow o Microsoft Dynamics spariranno. Penso invece che il loro ruolo cambierà.

Penso a [Claude Code](https://www.anthropic.com/research/claude-code-best-practices) e [GitHub Copilot](https://docs.github.com/en/copilot/get-started/about-github-copilot): già oggi possiamo delegare parte del lavoro di sviluppo a un agente. Copilot, per esempio, può analizzare un repository, modificarlo e preparare una pull request da revisionare. È l'inizio di un software più economico da adattare, non la fine dell'ingegneria, dei test o della manutenzione.

Oggi spesso adattiamo il lavoro alle schermate del software che abbiamo acquistato. Io immagino il contrario: il software che si adatta ai processi, costruito o modificato rapidamente a seconda delle necessità. Un'applicazione interna potrebbe vivere pochi mesi, essere sostituita da un'altra e lasciare intatti dati, API, identità e logiche di business.

La generazione di codice sta abbassando il costo di creare interfacce e automazioni. Gartner [stimava nel 2025](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027) che il 33% delle applicazioni software enterprise avrebbe incluso AI agentica entro il 2028, contro meno dell'1% nel 2024. Anche questa è una **proiezione**, non una misura della vita utile delle applicazioni.

La mia scommessa non è che tutto verrà riscritto continuamente. Quello sarebbe un incubo di manutenzione. È che **la parte visibile del software diventerà più economica da cambiare**, mentre contratti di dati, autorizzazioni, integrazioni, controlli e infrastrutture dovranno essere molto più stabili. In aziende regolamentate la velocità sarà comunque limitata da audit, sicurezza e affidabilità.

## 2028–2033 Smartphone, tablet e PC convergeranno

Penso che per gran parte delle attività quotidiane lo smartphone diventerà l'unico computer necessario: in tasca sarà un telefono, aperto sarà un tablet, collegato a uno schermo offrirà un desktop completo.

[Samsung DeX](https://www.samsung.com/it/apps/samsung-dex/) è uno degli esempi più vicini alla mia idea: su smartphone Galaxy compatibili, anche pieghevoli, permette di usare un ambiente desktop su uno schermo esterno. Un pieghevole con queste capacità assomiglia già a un computer da tenere in tasca. Il vero ostacolo rimane l'ecosistema software, non soltanto la potenza del dispositivo.

Non è fantascienza. Nel marzo 2026 [Google ha documentato](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays) le esperienze desktop di Android su monitor esterni, sviluppate anche con Samsung. Il problema non è più dimostrare che un telefono possa mostrare finestre su un monitor: è rendere l'esperienza sufficientemente solida per sostituire un portatile.

Rimangono ostacoli concreti: software professionale, periferiche, autonomia, dissipazione, multitasking pesante, architetture e compatibilità. Chi sviluppa videogiochi, fa rendering o lavora su ambienti particolari potrebbe continuare ad avere bisogno di workstation dedicate.

Non credo alla morte del PC. **Credo che il PC, per molte persone, diventerà una modalità di utilizzo anziché un dispositivo separato.**

## 2030–2037 Più protocolli condivisi, meno ecosistemi chiusi

Questa è probabilmente la previsione a cui tengo di più.

Il [Model Context Protocol (MCP)](https://www.anthropic.com/news/model-context-protocol), introdotto da Anthropic, è uno degli standard che trovo più interessanti. Immagino un servizio documentale o uno strumento GitHub esposto da un server MCP: più assistenti compatibili possono scoprire e utilizzare le operazioni che espone, nei limiti dei permessi concessi. [ChatGPT supporta app basate su MCP](https://help.openai.com/en/articles/12584461-developer-mode-and-full-mcp-connectors-in-chatgpt-beta) e Claude Code può collegarsi a server MCP. Questo non rende automaticamente portabili modelli, memoria e identità: standardizza il collegamento a strumenti e contesto. Vorrei che lo stesso principio si estendesse a molti altri servizi.

Io dovrei poter usare un'app di messaggistica, un'altra persona un client diverso, e dovremmo comunque comunicare. Come accade con le email. Lo stesso principio potrebbe estendersi a social network, contenuti, servizi digitali e AI.

I pezzi esistono già. [ActivityPub](https://www.w3.org/TR/activitypub/) è una raccomandazione W3C dal 2018 per il social networking decentralizzato. Il [Digital Markets Act europeo](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en) impone specifiche forme di interoperabilità ai grandi servizi di messaggistica designati: i primi servizi terzi interoperabili con WhatsApp sono arrivati nel novembre 2025. Nel mondo degli agenti, standard come MCP affrontano il problema di collegare strumenti e dati.

Sono segnali, non una garanzia. L'interoperabilità porta difficoltà di moderazione, antispam, crittografia, identità e incentivi economici. Le piattaforme hanno ottime ragioni commerciali per trattenere i propri utenti.

Credo però che il valore si sposterà gradualmente **dall'obbligo di usare la stessa app alla libertà di scegliere l'interfaccia che preferisci**.

## 2030–2038 Interfacce generate per quello che serve, quando serve

Oggi usiamo prodotti con schermate progettate in anticipo, più o meno uguali per tutti. Domani potremmo chiedere un risultato e ottenere un'interfaccia costruita attorno a quel compito: solo i dati utili, le azioni necessarie e i controlli adatti al dispositivo.

Le [app di ChatGPT](https://help.openai.com/en/articles/11487775-connectors-in-chatgpt) sono già un esempio concreto di servizi esterni portati dentro una conversazione. Il passo successivo che immagino è un'interfaccia composta al momento da componenti sicuri, in funzione del compito, anziché limitarsi a incorporare un'applicazione preesistente. È qui che progetti come A2UI diventano particolarmente interessanti.

Anche qui non partiamo da zero. [Google Research ha sperimentato interfacce generate dinamicamente](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), mentre il progetto aperto [A2UI](https://developers.googleblog.com/introducing-a2ui-an-open-project-for-agent-driven-interfaces/) definisce un modo per descrivere interfacce prodotte dagli agenti e renderizzarle su framework differenti.

Immagino di aprire la mia AI e chiederle di confrontare le spese di casa: invece di farmi navigare tra tre app, genera un grafico, una tabella e due azioni. Quando cambio obiettivo, cambia anche la schermata.

Questo non significa che ogni pixel debba essere improvvisato. Per operazioni delicate serviranno componenti affidabili, accessibilità, coerenza e conferme esplicite. **L'interfaccia può essere dinamica senza diventare imprevedibile.**

## 2030–2040 Sensori praticamente ovunque

Occhiali, cuffie, anelli, vestiti e oggetti della casa raccoglieranno sempre più informazioni sull'ambiente e sulle nostre attività. La mia AI personale potrebbe combinare alcuni di quei segnali per avere un contesto molto più utile rispetto a un semplice messaggio scritto.

Un esempio concreto è [Home Assistant con Matter](https://www.home-assistant.io/integrations/matter): dispositivi compatibili possono comunicare localmente senza dover passare sempre dal cloud del produttore. Se penso anche a dispositivi indossabili come i [Ray-Ban Meta](https://www.meta.com/ai-glasses/), immagino un assistente che sappia sfruttare il contesto della casa e dell'ambiente circostante. La loro integrazione sicura in un'unica AI personale è una mia ipotesi, non una funzionalità che oggi offrono insieme.

Un dato aiuta a capire la scala. Secondo [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/), nel mondo c'erano **18,5 miliardi** di dispositivi IoT connessi nel 2024; la stima per il 2025 era **21,1 miliardi** e la previsione per il 2030 arriva a **39 miliardi**.

![Grafico: connessioni IoT mondiali, 18,5 miliardi nel 2024, 21,1 stimati nel 2025, 39 previsti nel 2030.](/blog/charts/iot-devices.svg)

*Grafico 2 — [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/). Il totale comprende anche industria, infrastrutture e altri oggetti connessi: non sono tutti dispositivi indossabili. Il 2030 è una previsione.*

I sensori, però, non producono automaticamente informazioni utili. Dati imprecisi o fuori contesto possono essere peggiori dell'assenza di dati. E se un sistema conosce abitudini, posizione e segnali personali, privacy e sicurezza non possono essere opzioni secondarie.

Il futuro che vorrei è fatto di **controllo granulare, elaborazione locale quando possibile e condivisione minima dei dati**. Non di sorveglianza permanente mascherata da comodità.

## 2032–2040 Un'identità digitale che ci segue tra le piattaforme

Oggi account, preferenze, file e cronologie sono distribuiti tra servizi che spesso non comunicano. Immagino una maggiore separazione tra **chi sono**, **quali dati possiedo** e **quale applicazione uso per accedervi**.

L'Europa si sta muovendo in questa direzione. Il [wallet europeo di identità digitale](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en) dovrà essere reso disponibile dagli Stati membri entro la fine del 2026, con strumenti per identificarsi e condividere specifiche attestazioni. Il [Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act), applicabile dal settembre 2025, introduce anche misure per l'accesso ai dati di oggetti connessi e per cambiare fornitori cloud.

Attenzione, però: **un wallet d'identità non equivale a un profilo universale con tutte le nostre conversazioni, preferenze e memorie AI**. Quella portabilità più ampia è ancora da costruire, sia tecnicamente sia dal punto di vista giuridico.

La mia previsione è che sposteremo gradualmente il centro di gravità dall'account della piattaforma al controllo dell'utente. Non basterà un solo standard, e non succederà ovunque alla stessa velocità.

## 2035–2045 La realtà aumentata diventerà davvero spaziale

Non credo che il punto d'arrivo dell'AR sia avere notifiche e finestre 2D davanti agli occhi. Mi aspetto oggetti digitali che riconoscano la geometria del mondo, mantengano la loro posizione e reagiscano agli spazi reali.

[Apple Vision Pro e visionOS](https://developer.apple.com/visionos/) mostrano già come possano funzionare applicazioni spaziali. [Meta Orion](https://www.meta.com/design-at-meta/blog/three-design-leaders-share-their-vision-for-wearables-meta/), invece, è stato presentato come prototipo di occhiali AR. Sono approcci differenti: un visore disponibile come prodotto e un prototipo che esplora il formato degli occhiali. Nessuno dei due dimostra che l'AR spaziale, leggera e utilizzabile tutto il giorno, sia un problema risolto.

Immagino istruzioni ancorate a un oggetto mentre lo riparo, informazioni su ciò che sto osservando, o superfici virtuali collocate dove mi servono. È una differenza enorme rispetto a proiettare un piccolo schermo nel campo visivo.

La diffusione dei dispositivi non va confusa con la maturità di questa visione. [IDC](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/) prevede per la categoria degli occhiali con display circa **3 milioni di unità spedite nel 2026** e **12,2 milioni nel 2030**. Ma quella categoria comprende anche dispositivi che mostrano immagini davanti agli occhi senza offrire necessariamente AR spaziale completa.

I veri limiti saranno ancora autonomia, calore, peso, campo visivo, tracciamento preciso, qualità ottica e accettazione sociale. Per questo colloco la maturità della *vera* AR più lontano.

## 2040–2050+ Gli occhiali potrebbero sostituire lo smartphone

È la mia previsione più lontana e anche quella su cui ho più dubbi riguardo alle date.

Distinguerei i [Ray-Ban Meta](https://www.meta.com/ai-glasses/) dai [Meta Ray-Ban Display](https://www.meta.com/ai-glasses/meta-ray-ban-display/). I primi mostrano il potenziale dell'AI vocale, della fotocamera e delle interazioni a mani libere senza display; i secondi aggiungono un display nella lente e il controllo attraverso Meta Neural Band. Meta ha anche [documentato esperienze web per gli occhiali con display](https://developers.meta.com/blog/meta-connect-recap-ai-glasses/). Per me sono i primi passi verso una nuova interfaccia, non ancora un sostituto completo dello smartphone. Inoltre, mostrare informazioni nella lente non è la stessa cosa che avere oggetti digitali ancorati nello spazio reale.

Penso che a un certo punto gli occhiali potrebbero assorbire molte funzioni del telefono: comunicazione, navigazione, informazioni contestuali, traduzione, fotografie, interazione con la propria AI. Non avremmo più bisogno di tirare fuori uno schermo per ogni azione.

Il mercato si sta muovendo, anche se siamo lontani da una sostituzione. IDC prevede che gli smart glasses **senza display** passino da **13,6 milioni di spedizioni annue nel 2026** a **27,3 milioni nel 2030**. Gli occhiali con display costituiscono un mercato distinto e più piccolo.

![Grafico: spedizioni previste di occhiali intelligenti con e senza display nel 2026 e nel 2030.](/blog/charts/smart-glasses.svg)

*Grafico 3 — [IDC, 2026](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/). Si confrontano spedizioni annue previste, non utenti attivi né quota di sostituzione degli smartphone.*

I numeri dimostrano interesse e crescita prevista, **non dimostrano che lo smartphone scomparirà**. Per leggere a lungo, lavorare su contenuti complessi o usare un display privato potremmo continuare a preferire altri dispositivi.

La mia scommessa è più precisa: gli occhiali potrebbero diventare l'interfaccia primaria per gran parte delle interazioni brevi e contestuali. Il telefono potrebbe restare, ma passare in secondo piano.

## Il vero cambiamento: dal software al sistema che collega tutto

Se metto insieme queste nove previsioni, vedo un filo comune:

- **Le AI personali** potrebbero diventare il punto da cui partono molte attività.
- **Dati, identità e protocolli** dovrebbero rimanere utilizzabili anche cambiando software.
- **Le interfacce** potrebbero adattarsi al compito invece di imporre un percorso fisso.
- **I dispositivi** potrebbero diventare modi diversi di accedere allo stesso ambiente digitale.

Non credo a un futuro senza applicazioni. Credo a un futuro in cui le applicazioni saranno **meno centrali e più sostituibili**.

Anzi, il risultato potrebbe essere anche l'opposto di quello che spero: AI personali chiuse, dati ancora più concentrati, sensori invasivi e nuove dipendenze da pochi fornitori. La tecnologia da sola non garantisce apertura. Servono standard, regolamentazione, prodotti progettati bene e persone che chiedano di mantenere il controllo.

Per questo, tra tutte le previsioni, la più importante per me non riguarda un dispositivo o una data.

**Il vero valore non sarà soltanto nel software che utilizziamo, ma nella capacità dei sistemi di comunicare tra loro senza portarci via la libertà di scegliere.**

---

### Fonti e metodo

Ho usato dati pubblicati entro il **24 novembre 2026**. I numeri dei grafici provengono da fonti indicate accanto a ogni visualizzazione. Dove è presente una stima o una previsione di mercato l'ho etichettata come tale; le finestre temporali del post sono invece **ipotesi personali**, non proiezioni statistiche.

Fonti principali: [Stanford AI Index 2026](https://hai.stanford.edu/ai-index/2026-ai-index-report), [McKinsey — State of AI 2025](https://www.mckinsey.com/capabilities/quantumblack/our-insights/the-state-of-ai-2025), [Gartner — Agentic AI](https://www.gartner.com/en/newsroom/press-releases/2025-06-25-gartner-predicts-over-40-percent-of-agentic-ai-projects-will-be-canceled-by-end-of-2027), [IoT Analytics](https://iot-analytics.com/Number-Connected-Iot-Devices/), [IDC — Smart Glasses](https://www.idc.com/resource-center/blog/smart-glasses-surge-the-xr-market-is-rewriting-its-own-rules/), [Android Developers](https://developer.android.com/blog/posts/android-devices-extend-seamlessly-to-connected-displays), [W3C ActivityPub](https://www.w3.org/TR/activitypub/), [Linux Foundation — Agentic AI Foundation](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-Agentic-ai-foundation), [Google Research — Generative UI](https://research.google/blog/generative-ui-a-rich-custom-visual-interactive-user-experience-for-any-prompt/), [Commissione europea — Digital Identity](https://commission.europa.eu/topics/digital-economy-and-society/european-digital-identity_en), [Commissione europea — Data Act](https://digital-strategy.ec.europa.eu/en/policies/data-act) e [Digital Markets Act — Interoperability](https://digital-markets-act.ec.europa.eu/businesses-portal/messaging-interoperability_en).
