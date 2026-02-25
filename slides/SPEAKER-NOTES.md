# Speaker Notes — 10 Tips for Becoming a Top 1% AI User

Total talk time: 40-45 minutes + Q&A
24 slides. Average ~1.5 min per slide, with demo taking ~5 min.

---

## Slide 1: Title (2 min)

"Welcome everyone. I'm Manuel Corpas, Senior Lecturer at the University of Westminster and Turing Fellow. Tonight I want to share 10 practical techniques that separate casual AI users from the top 1% who are shipping research faster than ever."

"This is about AI fluency: the habits, tools, and workflows that compound your productivity. We'll go from quick wins you can adopt tomorrow morning, to building your own open-source bioinformatics skills. And I've got a live demo of something I built this week."

---

## Slide 2: The Gap (3 min)

"Here's the gap I see every day. Most researchers use AI like a search engine. They open ChatGPT, ask a question, paste the answer somewhere, and move on. No memory. No automation. No compounding."

"The top 1% have moved past that. They have agents that work while they sleep. Persistent memory that spans months. Automated pipelines that process papers, emails, and data without human intervention."

"The good news: you can cross this gap in a weekend. Let me show you how."

---

## Slide 3: Tip 1 (2 min)

"Tip 1 is the lowest-hanging fruit. Stop going to a browser to talk to AI. Put it inside your development environment."

"I use Claude Code in my terminal. It reads my entire codebase, understands the file structure, edits files directly, runs tests, and commits to git. It's not a chatbot. It's a pair programmer that knows your whole project."

"If you do one thing after tonight, install Claude Code."

---

## Slide 4: Tip 2 (2 min)

"Tip 2: build a prompt library. Every time you start a new Claude session, it reads a CLAUDE.md file at your project root. This is your AI onboarding document."

"Mine contains Python standards, architecture rules, path patterns, common commands. I never have to repeat instructions. The AI already knows how my project works."

"Treat your CLAUDE.md like you'd treat onboarding docs for a new team member. The better the docs, the better the AI performs."

---

## Slide 5: Tip 3 (2 min)

"Tip 3: memory. The biggest limitation of AI is that every conversation starts from zero. You solve this with a RAG pipeline."

"I embed my research notes, emails, publications, meeting transcripts, and session histories into ChromaDB. Then I can query semantically: 'What did I decide about the Wellcome proposal?' and get a real answer with source citations."

"I built my own RAG system for this. It is the foundation everything else builds on. The key insight: once your AI has memory, it stops being a chatbot and becomes a research partner."

---

## Slide 6: Tip 4 (1.5 min)

"Tip 4: use voice. Your brain outputs ideas faster as speech than typing. I run Whisper locally on Apple Silicon. Record a 5-minute voice memo, transcribe it, feed it to Claude for structuring."

"This is how I draft paper sections, plan projects, capture ideas on walks. Dictation is faster and often more natural than typed prompts."

---

## Slide 7: Tip 5 (2 min)

"Tip 5: automate your daily intelligence gathering. I have scheduled jobs that run every morning at 6:30am."

"arXiv ranking: my bot scrapes yesterday's preprints, ranks them by relevance to my research, and sends me the top 5 before breakfast."

"Podcast extraction: it transcribes and summarises AI podcasts I follow, extracting the key insights."

"Email triage: it classifies my inbox by urgency and drafts reply suggestions."

"By the time I open my laptop, the tedious work is already done."

---

## Slide 8: Tip 6 (2 min)

"Tip 6: persistent agents. I have two. RoboTerri runs on Telegram with 15 commands and 13 tools. It handles daily operations: paper summaries, podcast publishing, email drafts, writing assistance."

"RoboIsaac runs on WhatsApp as an analytical critique partner. He's modelled on Isaac Newton: rigorous, first-principles thinking, Socratic questioning. When I have an idea, I run it past Isaac first."

"Both share the same memory bridge. They never forget a conversation."

---

## Slide 9: Tip 7 (1.5 min)

"Tip 7: compounding. The top 1% don't just produce more output. They produce more output per insight."

"One research finding becomes: a paper paragraph, a blog post, a tweet thread, a podcast episode, a LinkedIn article. The AI handles the reformatting. I handle the thinking."

"This is how you go from publishing one paper a year to shipping continuously."

[PAUSE] "Now, tips 1 to 7 are about using AI. Tips 8 to 10 are about building with it. This is where it gets interesting for bioinformaticians."

---

## Slide 10: Transition (30 sec)

Just let this land. Brief pause.

"This is where it gets interesting for bioinformaticians."

---

## Slide 11: Tip 8 (3 min)

"Tip 8: contribute to open-source AI. Here's the problem. General purpose AI agents are powerful, but they're blind to biology."

"Three specific problems:"

"First: privacy. Genomic data is sensitive. You cannot send your patient VCFs to a cloud API. We need local-first execution."

"Second: reproducibility. Biology demands audit trails. Every analysis step must be logged, versioned, and exportable as a reproducible pipeline."

"Third: domain knowledge. A generic agent doesn't know that a VCF file needs ancestry-aware annotation, or that single-cell data needs doublet removal before clustering."

[PAUSE] "So I'm announcing tonight: ClawBio. The first bioinformatics-native AI agent skill library. Open source. Local first. Privacy focused."

---

## Slide 12: Architecture (2 min)

"Here's how it works. You describe what you want in natural language. The Bio Orchestrator detects your file type, routes to the right specialist skill, runs the analysis, and produces a markdown report with figures, tables, and a reproducibility bundle."

"Each skill wraps proven bioinformatics tools. Biopython, SAMtools, Scanpy, AlphaFold. The AI orchestrates; the tools compute."

---

## Slide 13: Skill Catalogue (1.5 min)

"Eight skills planned. Two are at MVP: the Equity Scorer and the Bio Orchestrator. Six more are on the roadmap."

"Each skill is just a SKILL.md file plus Python scripts. Modular. Composable. You can use one alone or chain them through the orchestrator."

"Let me show you what the Equity Scorer does."

---

## Slide 14: Tip 9 + Demo Intro (switch to terminal)

"Tip 9: build modular skills. Let me run one live."

[SWITCH TO TERMINAL — see DEMO-SCRIPT.md]

---

## Slides 15-19: Demo Results (5 min total)

Walk through terminal output, PCA, FST heatmap, ancestry chart, HEIM score.
See DEMO-SCRIPT.md for detailed talking points.

---

## Slide 20: Why This Matters (2 min)

"Why does this matter? 86% of GWAS participants are European. That means polygenic risk scores, drug targets, and clinical guidelines are biased towards one population."

"The HEIM Index gives researchers a single number to quantify this problem. Score your dataset. Report it alongside your demographics. Track it over time."

"I have a paper in review on this. But the tool is open source, and you can run it tonight."

---

## Slide 21: Tip 10 + CTA (2 min)

"Tip 10: create infrastructure that ships your research. The Equity Scorer took 2 days to build. Real population genetics, not a toy. Anyone in this room can build the next skill."

"I need skills from this community. GWAS pipeline? Metagenomics classifier? Clinical variant reporter? Pathway enricher? The template is there. The orchestrator routes to your skill automatically."

"Who wants to build one?" [Look around, make eye contact]

---

## Slide 22: Get Started (1.5 min)

"Here's the link. The repo goes live tonight. Star it, clone it, run the Equity Scorer on your own data."

"If you want to contribute a skill, the template and contributing guide are in the repo. Open a PR. I'll review it personally. It's MIT licensed, local-first, and designed so that anyone in this room can add to it."

---

## Slide 23: Recap (1 min)

Quick recap of all 10 tips. Don't linger; just hit the main line:

"Tips 1 to 7: AI fluency, the productivity habits that compound. Tips 8 to 10: building with AI, turning your domain expertise into open-source tools. The throughline: the top 1% build systems that compound. Everyone else uses tools one at a time."

---

## Slide 24: Thank You

"Thank you. I'll be at the pub after. Come talk to me if you want to build a skill, have ideas for the library, or just want to chat about AI in bioinformatics."

"Questions?"

---

## Timing Summary

| Section | Slides | Time |
|---------|--------|------|
| Title + Gap | 1-2 | 5 min |
| Tips 1-7 | 3-9 | 13 min |
| Transition | 10 | 0.5 min |
| Tips 8 (Announcement) | 11-13 | 6.5 min |
| Tip 9 (Live Demo) | 14-19 | 8 min |
| Why It Matters | 20 | 2 min |
| Tip 10 + CTA | 21-22 | 3.5 min |
| Recap + Thank You | 23-24 | 2 min |
| **Total** | **24** | **~40 min** |
| Q&A | — | 5-15 min |
