# SEO operativo — donatomolino.it

Questo documento contiene i passaggi che richiedono il tuo accesso agli account. Le modifiche al codice sono già nella PR.

## 1. Dopo il merge e il deploy

Verifica da una finestra anonima gli URL:
- https://donatomolino.it/
- https://donatomolino.it/blog/it/
- https://donatomolino.it/blog/en/
- https://donatomolino.it/blog/it/il-futuro-del-tech-2027-2050/
- https://donatomolino.it/blog/en/the-future-of-tech-2027-2050/
- https://donatomolino.it/sitemap.xml
- https://donatomolino.it/robots.txt
- https://donatomolino.it/llms.txt
- https://donatomolino.it/blog/rss.xml

La sitemap include automaticamente nuovi articoli pubblicati ad ogni build.

## 2. Google Search Console
1. Vai a https://search.google.com/search-console/
2. Aggiungi la proprietà di tipo Dominio per donatomolino.it.
3. Verifica la proprietà nel DNS: copia il valore TXT esatto fornito da Google nelle impostazioni del tuo provider DNS (ad esempio Cloudflare se gestisce il dominio).
4. Nella sezione Sitemap invia https://donatomolino.it/sitemap.xml
5. Usa Controllo URL sulle due versioni dell'articolo; se appropriato, richiedine l'indicizzazione.
6. Verifica nel tempo rapporti su indicizzazione, clic/impressioni e Core Web Vitals. L'indicizzazione non è garantita.

## 3. Bing Webmaster Tools
1. Apri https://www.bing.com/webmasters/
2. Importa il dominio verificato da Google Search Console, oppure verifica manualmente.
3. Controlla sitemap e URL indicizzati.

## 4. Verifica sicurezza e accesso bot
- robots.txt permette OAI-SearchBot, ChatGPT-User, Claude-SearchBot e Claude-User, oltre ai crawler generici.
- Se usi Cloudflare o un altro WAF, controlla AI Crawl Control, protezioni bot, firewall, blocchi per IP o challenge. Un permesso robots.txt non supera una protezione sul server.
- Crawling per ricerca e crawling per addestramento sono finalità distinte; puoi regolarle separatamente.

## 5. Contatore lettori (opzionale)
Nel tuo account Cloudflare puoi abilitare Web Analytics (o l'opzione one-click su Pages, se presente).
https://dash.cloudflare.com/
Le visualizzazioni e i visitatori unici sono stime di traffico e non letture complete certificate.

## Precisazioni
- llms.txt è un riferimento supplementare: Google non lo richiede e non lo usa per AI Overviews.
- Sitemap, JSON-LD e metadati non garantiscono ranking, indicizzazione o citazioni nelle risposte AI.
- Lavorare sulla qualità dei contenuti originali e acquisire link editoriali pertinenti conta più delle scorciatoie SEO.
