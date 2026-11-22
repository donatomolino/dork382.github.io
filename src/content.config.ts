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

const blog = defineCollection({
  loader: glob({
    base: "./src/content/blog",
    pattern: "**/*.{md,mdx}",
  }),
  schema: z.object({
    title: z.string().min(1),
    description: z.string().min(1),
    pubDate: z.coerce.date(),
    updatedDate: z.coerce.date().optional(),
    draft: z.boolean().default(false),
    tags: z.array(z.string()).default([]),
  }),
});

export const collections = {
  showcase,
  blog,
};
