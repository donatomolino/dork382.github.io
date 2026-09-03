import icon from "astro-icon";
import react from "@astrojs/react";
import { fileURLToPath } from "node:url";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://donatomolino.it",
  vite: {
    resolve: { alias: { "~": fileURLToPath(new URL("./src", import.meta.url)) } },
  },
  devToolbar: {
    enabled: false,
  },
  integrations: [icon(), react()],
});
