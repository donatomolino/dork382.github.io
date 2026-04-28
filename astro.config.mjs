import icon from "astro-icon";
import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://donatomolino.it",
  devToolbar: {
    enabled: false,
  },
  integrations: [icon()],
});
