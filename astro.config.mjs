import icon from "astro-icon";
import react from "@astrojs/react";
import node from "@astrojs/node";
import { fileURLToPath } from "node:url";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://donatomolino.it",
  output: "server",
  adapter: node({ mode: "standalone" }),
  vite: {
    resolve: { alias: { "~": fileURLToPath(new URL("./src", import.meta.url)) } },
  },
  devToolbar: {
    enabled: false,
  },
  integrations: [icon(), react()],
});
