CREATE TABLE IF NOT EXISTS comments (
  id TEXT PRIMARY KEY,
  thread TEXT NOT NULL,
  parent_id TEXT REFERENCES comments(id),
  name TEXT NOT NULL CHECK(length(name) BETWEEN 1 AND 60),
  body TEXT NOT NULL CHECK(length(body) BETWEEN 1 AND 2000),
  lang TEXT NOT NULL CHECK(lang IN ('it','en')),
  status TEXT NOT NULL DEFAULT 'pending' CHECK(status IN ('pending','approved','rejected')),
  created_at INTEGER NOT NULL,
  day TEXT NOT NULL,
  rate_key TEXT
);
CREATE INDEX IF NOT EXISTS comments_thread_status ON comments(thread,status,parent_id,created_at);
CREATE INDEX IF NOT EXISTS comments_day_rate ON comments(day,rate_key);
