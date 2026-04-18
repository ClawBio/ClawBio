# Pointing a Custom Domain at an OpenClaw Webchat

An **OpenClaw gateway** is the self-hosted conversational layer that sits between a browser chat client and your local ClawBio runtime. The browser sends a natural-language request to the gateway, the gateway asks an external LLM to interpret it, and that LLM translates the request into local ClawBio CLI or Python calls.

This is the plain-language alternative to Telegram or Discord adapters. You might choose it if you want browser-based chat without relying on messenger platforms, if you want tighter operational control over authentication and deployment, or if you want the strongest privacy story while keeping ClawBio's biological processing local.

## What you need

- A domain you own (e.g. `yourdomain.com`)
- An OpenClaw gateway running with Tailscale Funnel enabled
- Your Tailscale Funnel URL (looks like `https://<machine-name>.<tailnet>.ts.net`)

## Option A: Cloudflare Proxy (Recommended)

Cloudflare (free tier) handles TLS and proxies requests to your Tailscale Funnel address. The user sees your clean domain; the `.ts.net` address stays hidden.

### 1. Add your domain to Cloudflare

1. Sign up at https://dash.cloudflare.com (free)
2. Click **Add a site** and enter your domain
3. Cloudflare will scan your existing DNS records — review and keep them
4. Cloudflare will give you two nameservers (e.g. `anna.ns.cloudflare.com`, `bob.ns.cloudflare.com`)

### 2. Update nameservers at your registrar

Go to your domain registrar (GoDaddy, Namecheap, Google Domains, Porkbun, etc.):

1. Find **DNS** or **Nameserver** settings
2. Switch from default nameservers to **Custom**
3. Enter the two Cloudflare nameservers
4. Save — propagation usually takes 10 minutes to 1 hour

### 3. Add a CNAME record in Cloudflare

In Cloudflare DNS settings, add a record:

| Field | Value |
|---|---|
| Type | `CNAME` |
| Name | `chat` (or `@` for root domain) |
| Target | Your Tailscale Funnel hostname, e.g. `my-machine.tail1234.ts.net` |
| Proxy status | **Proxied** (orange cloud ON) |

### 4. Set SSL/TLS mode

In Cloudflare: **SSL/TLS** → set encryption mode to **Full**.

This tells Cloudflare to connect to your Tailscale Funnel over HTTPS (Tailscale provides a valid cert on its end).

### 5. Update OpenClaw allowed origins

Edit your OpenClaw config (`~/.openclaw/openclaw.json`) and add your new domain to the allowed origins:

```json
"controlUi": {
  "enabled": true,
  "allowedOrigins": [
    "https://<your-tailscale-hostname>.ts.net",
    "https://chat.yourdomain.com"
  ]
}
```

Then restart the gateway:

```bash
# macOS (launchd)
launchctl kickstart -k gui/$(id -u)/ai.openclaw.gateway

# Linux (systemd)
sudo systemctl restart openclaw-gateway

# Or just kill and relaunch manually
```

### 6. Test it

Visit `https://chat.yourdomain.com` — you should see the OpenClaw Control UI or your webchat page.

---

## Option B: Registrar Redirect (No Cloudflare)

If you don't want to move DNS, most registrars support subdomain forwarding:

1. Go to your registrar's DNS or Forwarding settings
2. Add a subdomain forward:
   - **Subdomain**: `chat`
   - **Forward to**: `https://<your-tailscale-hostname>.ts.net/`
   - **Type**: Temporary (302) or Permanent (301)

This works but the URL bar will change to the `.ts.net` address after redirect.

---

## Option C: Static-hosted webchat with custom domain

Host the webchat HTML on a static hosting service (GitHub Pages, Vercel, Netlify, Cloudflare Pages) and have it connect to your gateway over WebSocket.

### 1. Deploy the webchat HTML

Put `chat/index.html` in a GitHub repo, Vercel project, or similar. Set the `GATEWAY_URL` in the HTML to your Tailscale Funnel WebSocket address:

```js
const GATEWAY_URL = 'wss://<your-tailscale-hostname>.ts.net';
```

### 2. Set up the custom domain

Follow your hosting provider's docs to point `chat.yourdomain.com` at the deployment:

- **GitHub Pages**: Add a `CNAME` file containing `chat.yourdomain.com`, then add a CNAME DNS record pointing `chat` to `<username>.github.io`
- **Vercel**: Add the domain in project settings, then add a CNAME DNS record pointing `chat` to `cname.vercel-dns.com`
- **Netlify**: Add the domain in site settings, then add a CNAME DNS record

### 3. Update OpenClaw allowed origins

Same as Option A step 5 — add `https://chat.yourdomain.com` to `controlUi.allowedOrigins`.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| "unauthorized: gateway password missing" | Enter the gateway password in the Control UI settings panel, or pass it as `?p=yourpassword` in the URL |
| "This device needs pairing approval" | On the gateway host, run `openclaw devices list` then `openclaw devices approve <requestId>` |
| WebSocket connection fails | Check that Tailscale Funnel is enabled: `tailscale funnel status`. Ensure `gateway.tailscale.mode` is `"funnel"` in openclaw.json |
| CORS / origin blocked | Add your domain to `controlUi.allowedOrigins` in openclaw.json and restart the gateway |
| SSL errors with Cloudflare | Make sure SSL/TLS mode is set to **Full** (not Flexible or Full Strict) |
