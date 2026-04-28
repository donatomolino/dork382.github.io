/// <reference path="../.astro/types.d.ts" />
/// <reference types="astro/client" />

interface ImportMetaEnv {
  readonly TELEGRAM_CHANNEL_USERNAME: string
}

interface ImportMeta {
  readonly env: ImportMetaEnv
}
