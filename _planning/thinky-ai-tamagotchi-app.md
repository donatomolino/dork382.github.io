# Thinky AI Tamagotchi App

Documento di ragionamento per trasformare Thinky da esperienza web interattiva a companion app mobile in stile tamagotchi con AI.

Questa cartella e' fuori da `src` e `public`, quindi non viene inclusa nel deploy Astro.

## Visione

Thinky diventa una piccola presenza quotidiana: non solo una mascotte che risponde, ma un compagno vivo che ha bisogni leggeri, memoria, umore, abitudini e una relazione con una persona, una coppia o un piccolo gruppo.

La versione attuale del sito e' un ottimo prototipo di identita': design, mood, micro-animazioni, chat sintetica, personalita'. L'app dovrebbe mantenere quella semplicita' e aggiungere un loop quotidiano.

Obiettivo emotivo: aprire l'app deve sembrare come controllare un piccolo essere digitale, non gestire una dashboard.

## Loop Core

Ogni giorno Thinky ha 3 bisogni principali:

- Food: va nutrito una volta al giorno.
- Water: deve bere una volta al giorno.
- Play: vuole giocare/interagire una volta al giorno.

Ogni bisogno influenza:

- mood visivo;
- frase brevissima;
- stato nel diario;
- streak giornaliero;
- qualita' della notifica successiva;
- relazione con gli adottanti.

Regola importante: niente punizioni pesanti. Thinky non deve morire o far sentire in colpa. Se viene trascurato, diventa piu' pensieroso, pigro, buffo o bisognoso, ma resta accogliente.

## Routine Giornaliera

Una giornata ideale:

1. Al mattino Thinky manda una notifica spontanea.
2. L'utente apre l'app.
3. Thinky mostra il suo stato: fame, sete, voglia di giocare, mood.
4. L'utente fa 1-3 azioni rapide.
5. Thinky risponde con una frase breve e un'animazione.
6. A fine giornata Thinky scrive una riga di diario.

Esempio:

- Notifica: "I found a tiny thought."
- Azione: l'utente gli da' acqua.
- Risposta: "Glug. Better."
- Diario: "Today I drank late, but I stayed cheerful."

## Diario Personale

Thinky dovrebbe avere un diario generato una volta al giorno.

Il diario non deve essere lungo. Deve sembrare scritto da lui:

- 1-3 frasi;
- tono tenero, sintetico, un po' curioso;
- riferimenti alle azioni reali della giornata;
- niente over-sharing artificiale;
- niente "AI assistant vibe".

Esempi:

- "Donato fed me a cookie. I glowed blue after that. Good day."
- "Nobody visited for a while. I counted grid lines and waited."
- "Two humans checked on me today. I felt shared."

Il diario puo' diventare il cuore della retention: dopo 30 giorni l'utente ha una piccola storia affettiva.

## Notifiche

Thinky dovrebbe scrivere una volta al giorno, non bombardare.

Tipi di notifica:

- Care reminder: "I may be tiny-hungry."
- Diary teaser: "I wrote one line today."
- Mood ping: "My glow changed."
- Group prompt: "Someone should give me water."
- Streak gentle nudge: "Day 6 is waiting."

Criticita':

- iOS richiede permesso esplicito e APNs.
- Android ha permessi moderni per notifiche.
- Le notifiche vanno schedulate lato backend per affidabilita'.
- Notifiche locali possono bastare per MVP, ma non funzionano bene per messaggi AI spontanei sincronizzati su coppie/gruppi.

Fonti tecniche: Expo Notifications gestisce token, notifiche locali e remote in React Native/Expo; Firebase Cloud Messaging e' lo standard cross-platform per push. FlutterFire/Firebase copre anche Flutter, ma iOS richiede comunque setup APNs e test su device fisico.

## Streak

La streak non dovrebbe essere solo "hai aperto l'app".

Streak valida se almeno una di queste e' vera:

- Thinky ha mangiato;
- Thinky ha bevuto;
- Thinky ha giocato;
- qualcuno del gruppo ha fatto una care action.

Tipi:

- Personal streak: giorni in cui tu hai interagito.
- Household streak: giorni in cui il gruppo si e' preso cura di Thinky.
- Bond streak: giorni consecutivi in cui Thinky ha scritto diario.

Idea: se si perde la streak, Thinky puo' conservare una "memory scar" morbida, tipo:

"We skipped day 12. I kept the spot warm."

## Adozione Singola, Coppia, Gruppo

Questa feature e' forte.

Modelli:

- Single adoption: una persona, un Thinky.
- Couple adoption: due persone, un Thinky condiviso.
- Group adoption: 3-6 persone, un Thinky condiviso.

Ogni Thinky ha:

- owner principale;
- members;
- ruoli leggeri;
- diario condiviso;
- stato condiviso;
- azioni giornaliere attribuite.

Esempio dati:

```txt
thinky_households
- id
- name
- created_by
- mode: single | couple | group
- current_mood
- hunger_state
- thirst_state
- play_state
- streak_count
- last_daily_reset_at

members
- user_id
- household_id
- role: owner | member
- display_name
- joined_at

care_events
- id
- household_id
- user_id
- type: feed | drink | play | chat
- metadata
- created_at

diary_entries
- id
- household_id
- date
- text
- mood
- generated_from_event_ids
```

Criticita' multiutente:

- conflitti: due persone danno da mangiare nello stesso momento;
- notifiche duplicate;
- privacy del diario;
- inviti e uscita dal gruppo;
- moderazione del testo se chat condivisa.

Soluzione MVP:

- un solo Thinky per household;
- invito tramite link;
- massimo 4 membri;
- azioni giornaliere idempotenti: se ha gia' mangiato, seconda azione diventa snack/playful extra.

## Login

Login con Google o Apple ha senso da subito se l'app e' mobile e multiutente.

Per iOS, Sign in with Apple diventa praticamente obbligatorio se offri altri login social.

Opzioni backend auth:

- Firebase Auth: rapido, solido, Google/Apple, integrazione con FCM.
- Supabase Auth: buono se vuoi Postgres, RLS, SQL e un backend piu' leggibile.
- Custom auth: non consigliato per MVP.

Preferenza pratica:

- Firebase Auth + Firestore per MVP veloce.
- Supabase se vuoi piu' controllo dati, SQL, relazioni, export e query complesse.

Per Thinky, io sceglierei Supabase se vuoi costruire un prodotto duraturo e leggibile. Sceglierei Firebase se vuoi arrivare prima a notifiche e auth mobile senza troppe frizioni.

## AI Design

L'AI non deve essere una chat generica. Deve essere un motore di micro-personalita'.

Output consigliato:

```json
{
  "text": "Tiny answer.",
  "mood": "curioso",
  "action": "cookie",
  "diary_hint": "Donato fed me today.",
  "care_signal": "fed"
}
```

Regole:

- risposta sempre brevissima;
- massimo 80 caratteri in bubble;
- diario massimo 220 caratteri;
- lingua dell'utente;
- niente spiegazioni lunghe;
- niente consigli medici/legali/finanziari;
- se non capisce, risponde comunque con tono Thinky.

Uso dell'AI:

- risposta chat;
- diario giornaliero;
- notifica giornaliera;
- riassunto memoria settimanale;
- reazione a care action.

Non usare AI per:

- calcolare fame/sete/streak;
- decidere permessi;
- gestire sicurezza;
- stato canonico del pet.

La logica di stato deve essere deterministica. L'AI deve colorare l'esperienza, non comandarla.

## Memoria

Memoria utile:

- nome degli utenti;
- routine preferita;
- eventi recenti;
- tono preferito;
- care pattern;
- inside jokes approvati;
- diario storico.

Memoria da evitare:

- dati sensibili;
- inferenze emotive troppo forti;
- contenuti privati non necessari;
- memoria infinita senza controllo.

Strategia:

- memoria breve: ultimi 10 eventi;
- memoria giornaliera: diario;
- memoria lunga: summary settimanale generato;
- possibilita' di cancellare memoria.

## Feature List

MVP:

- app mobile con Thinky animato;
- login Google/Apple;
- single household;
- feed/drink/play una volta al giorno;
- chat breve con Thinky;
- notifiche giornaliere;
- diary entry giornaliera;
- streak;
- basic profile.

MVP Plus:

- coppia/gruppo con invite link;
- azioni attribuite ai membri;
- diario condiviso;
- notifiche intelligenti per chi non ha ancora interagito;
- personalizzazione nome Thinky;
- alcuni oggetti sbloccabili.

V1:

- marketplace cosmetico leggero;
- mood history;
- recap settimanale;
- widget iOS/Android;
- lock screen/live activity semplice;
- streak freeze gentile;
- piu' minigiochi.

V2:

- voice notes;
- AR desk mascot;
- physical merch tie-in;
- desktop companion;
- multiple Thinky personalities;
- shared rituals per coppie.

## Minigiochi

Devono durare 10-30 secondi.

Idee:

- tap bubbles: Thinky guarda le bolle.
- feed timing: scegli il cibo giusto.
- glow match: abbina colore a mood.
- tiny puzzle: ricomponi una piccola idea.
- breathing game: Thinky si calma con te.
- doodle: disegni una forma e Thinky reagisce.

I minigiochi non devono diventare la parte principale. Sono rituali, non livelli.

## Economia e Retention

Non trasformarlo in un free-to-play aggressivo.

Possibili monetizzazioni:

- premium una tantum;
- abbonamento leggero per gruppi, piu' memoria, piu' diario;
- cosmetic packs;
- physical Thinky desk toy in futuro;
- B2C companion app gratuita con limite AI.

Evita:

- dark pattern;
- streak guilt;
- notifiche insistenti;
- paywall su bisogni base.

## Flutter vs React Native

React Native / Expo:

Pro:

- piu' vicino all'attuale stack React;
- riuso mentale di componenti, logica, TypeScript;
- Expo accelera prototipazione, build, notifiche, OTA updates;
- ottimo per MVP;
- piu' facile condividere parte della logica con web Thinky.

Contro:

- animazioni complesse richiedono cura con Reanimated/Skia;
- quando tocchi moduli nativi avanzati devi gestire EAS/native config;
- performance molto buona, ma serve disciplina.

Flutter:

Pro:

- UI e animazioni molto solide e coerenti;
- ottimo controllo grafico;
- performance fluida;
- buona scelta per un pet molto custom, animato, illustrativo.

Contro:

- stack diverso da quello attuale;
- meno riuso del codice React/TypeScript;
- integrazione AI/backend ok, ma il team deve amare Dart;
- piu' costo iniziale se parti da zero.

Scelta consigliata:

Per Thinky io partirei con React Native + Expo.

Motivo: hai gia' un Thinky React, hai gia' logica mood/action, e il prodotto deve validare retention, notifiche e diario piu' che grafica estrema. Expo ti fa arrivare prima a un'app testabile. Se poi la mascotte diventa molto piu' animata, puoi usare Skia/Reanimated o valutare Flutter dopo una validazione reale.

## Stack Suggerito

Opzione A, MVP veloce:

- React Native + Expo
- TypeScript
- Expo Notifications
- Firebase Auth
- Firestore
- Cloud Functions o Cloudflare Worker per AI
- OpenRouter per modello economico

Opzione B, prodotto piu' controllabile:

- React Native + Expo
- TypeScript
- Supabase Auth
- Supabase Postgres
- Edge Functions o Cloudflare Worker per AI
- FCM/APNs per notifiche

Io farei:

- Expo per app;
- Supabase per dati relazionali;
- Cloudflare Worker per AI e scheduled jobs;
- Firebase/FCM solo se serve un'integrazione push piu' diretta.

Nota: Expo Notifications puo' gestire push in modo pratico, ma per notifiche programmate server-side e produzione bisogna disegnare bene token, permessi e retry.

## Backend

Servizi necessari:

- auth;
- household management;
- pet state;
- care events;
- diary generation;
- push scheduling;
- AI proxy;
- rate limiting;
- moderation/safety;
- analytics minimale.

Job giornalieri:

- reset bisogni;
- generazione notifica;
- generazione diario;
- streak update;
- cleanup token push invalidi.

Possibile architettura:

```txt
Mobile App
  -> Supabase Auth
  -> Supabase DB
  -> Cloudflare Worker API
       -> OpenRouter
       -> Push provider
       -> scheduled daily jobs
```

## Criticita'

Notifiche:

- permessi non sempre concessi;
- su iOS serve device fisico per test affidabile;
- token possono scadere;
- fusi orari complicano "una volta al giorno".

AI:

- costi se cresce;
- risposte non JSON;
- finish_reason length;
- prompt injection;
- contenuti non desiderati in chat condivisa.

Soluzioni AI:

- schema rigido;
- max output piccolo ma sufficiente;
- fallback locale;
- retry una volta;
- logging anonimo degli errori;
- filtro di lunghezza e mood/action validi.

Multiutente:

- sincronizzazione realtime;
- azioni duplicate;
- gestione inviti;
- privacy e uscita dal gruppo.

Prodotto:

- se il loop e' troppo semplice, stanca;
- se e' troppo pesante, diventa un dovere;
- la magia sta nel diario e nelle micro-reazioni.

## Roadmap

Fase 0: Concept

- definire stati pet;
- definire 20 care reactions;
- definire tono AI;
- prototipo Figma o direttamente Expo.

Fase 1: Prototype locale

- Thinky animato in Expo;
- feed/drink/play locali;
- streak locale;
- diario finto;
- no login.

Fase 2: MVP backend

- login;
- household single;
- salvataggio stato;
- Worker AI;
- diario reale;
- notifiche base.

Fase 3: Closed beta

- 20-50 utenti;
- test retention 7 giorni;
- test notifiche;
- raccolta frasi belle/brutte;
- costo AI per utente.

Fase 4: Couple/group

- invite link;
- shared diary;
- attribution;
- notification routing.

Fase 5: Public launch

- polish;
- App Store / Play Store;
- privacy policy;
- terms;
- analytics;
- crash reporting.

## Metriche

Metriche sane:

- D1, D7, D30 retention;
- quanti utenti completano feed/drink/play;
- quante diary entry vengono lette;
- opt-in notifiche;
- risposte AI fallite;
- costo AI per active user;
- numero di household multiutente.

Metriche da non idolatrare:

- tempo infinito in app;
- notifiche aperte a tutti i costi;
- streak come pressione.

## Domande Aperte

- Thinky deve avere un nome personalizzabile o restare sempre Thinky?
- Deve poter "ammalarsi" o solo cambiare mood?
- Il diario e' privato o condiviso sempre?
- In coppia/gruppo, tutti vedono tutte le chat?
- Le notifiche devono sembrare scritte live o chiaramente programmate?
- Quanto deve costare l'AI per utente al mese?
- Ci saranno oggetti/collezionabili?
- L'app deve essere solo mobile o anche web login?

## Decisione Consigliata

Costruire un MVP React Native + Expo.

Prima milestone:

- single-user;
- no gruppo;
- feed/drink/play;
- 7-day streak;
- diario AI giornaliero;
- una notifica al giorno;
- login opzionale o anonimo.

Seconda milestone:

- login Google/Apple;
- sync cloud;
- coppia/gruppo.

Motivo: il rischio principale non e' tecnico, e' capire se le persone vogliono davvero tornare ogni giorno da Thinky. Prima validiamo il rituale, poi aggiungiamo la parte sociale.

## Fonti Tecniche Consultate

- Expo Notifications: https://docs.expo.dev/push-notifications/push-notifications-setup/
- Expo Notifications SDK: https://docs.expo.dev/versions/latest/sdk/notifications/
- Firebase Cloud Messaging: https://firebase.google.com/docs/cloud-messaging
- Firebase Authentication: https://firebase.google.com/docs/auth
- Firebase Flutter setup / APNs note: https://firebase.google.com/docs/flutter/setup
- FlutterFire FCM Apple integration: https://firebase.flutter.dev/docs/messaging/apple-integration/
