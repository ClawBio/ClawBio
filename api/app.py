"""ClawBio API — FastAPI application factory."""

from fastapi import FastAPI

from api.routers import jobs


def create_app() -> FastAPI:
    app = FastAPI(
        title="ClawBio API",
        description=(
            "REST API for ClawBio bioinformatics skills. "
            "Submit skill runs, track progress, and resume interrupted jobs."
        ),
        version="0.1.0",
    )
    app.include_router(jobs.router)

    @app.get("/health", tags=["meta"])
    async def health():
        return {"status": "ok"}

    return app


app = create_app()
