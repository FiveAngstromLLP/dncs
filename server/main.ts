import { Application } from "https://deno.land/x/oak/mod.ts";

const app = new Application();

app.use(async (ctx) => {
  if (ctx.request.url.pathname === "/api/generate") {
    if (ctx.request.method === "POST") {
      try {
        const body = await ctx.request.body().value;
        const config = JSON.parse(body);

        // Convert config to TOML format
        const tomlContent = Object.entries(config)
          .map(([key, value]) => `${key} = "${value}"`)
          .join("\n");

        // Generate a unique filename
        const filename = `config_${Date.now()}.toml`;
        const filePath = `${Deno.cwd()}/configs/${filename}`;

        // Ensure the configs directory exists
        await Deno.mkdir(`${Deno.cwd()}/configs`, { recursive: true });

        // Write the TOML file
        await Deno.writeTextFile(filePath, tomlContent);

        ctx.response.body = JSON.stringify({
          message: "Configuration saved successfully",
          file_path: filePath,
        });
        ctx.response.type = "application/json";
      } catch (error) {
        ctx.response.status = 400;
        ctx.response.body = JSON.stringify({ error: error.message });
        ctx.response.type = "application/json";
      }
    } else {
      ctx.response.status = 405;
      ctx.response.body = "Method Not Allowed";
    }
  } else {
    await ctx.send({
      root: `${Deno.cwd()}/static`,
      index: "index.html",
    });
  }
});

console.log("Server running on http://localhost:8000");
await app.listen({ port: 8000 });
