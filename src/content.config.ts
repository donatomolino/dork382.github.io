import { defineCollection } from "astro:content";
import { glob } from "astro/loaders";
import { z } from "astro/zod";

const isValidShowcaseUrl = (value: string) =>
  value.startsWith("/") || z.url().safeParse(value).success;

const showcase = defineCollection({
  loader: glob({
    base: "./src/content/showcase",
    pattern: "**/*.json",
  }),
  schema: ({ image }) =>
    z.object({
      title: z.string().min(1),
      image: image().optional(),
      url: z.string().min(1).refine(isValidShowcaseUrl, {
        error: "Expected an absolute URL or a root-relative path.",
      }),
      featured: z.number().min(1).optional(),
    }),
});

export const collections = {
  showcase,
};