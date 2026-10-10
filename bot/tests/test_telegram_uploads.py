"""Telegram upload handlers reject bad files with a reply or a silent return.

A NameError inside either handler is swallowed by its except block and reaches
the user as "Sorry, something went wrong -- NameError", so these run the real
handlers with a mocked Update rather than reading the source.
"""

import asyncio
import sys
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock

import pytest

pytest.importorskip("telegram")  # the bot dependency group; security tests stay stdlib-only

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bot/

import roboterri  # noqa: E402


def _update(*, document=None):
    update = MagicMock()
    update.effective_chat.id = 1
    update.message.photo = []
    update.message.document = document
    update.message.caption = None
    update.message.reply_text = AsyncMock()
    return update


def _document(file_name: str, mime_type: str):
    doc = MagicMock()
    doc.file_name = file_name
    doc.mime_type = mime_type
    doc.file_size = 10
    doc.get_file = AsyncMock(return_value=MagicMock(download_as_bytearray=AsyncMock(return_value=bytearray(b"x"))))
    return doc


@pytest.fixture
def context(monkeypatch):
    monkeypatch.setattr(roboterri, "_check_rate_limit", lambda update: True)
    ctx = MagicMock()
    ctx.bot.send_chat_action = AsyncMock()
    return ctx


def test_document_with_disallowed_extension_gets_the_unsupported_reply(context):
    update = _update(document=_document("payload.py", "text/x-python"))
    asyncio.run(roboterri.handle_document(update, context))
    reply = update.message.reply_text.await_args.args[0]
    assert reply.startswith("Unsupported file type (.py)")
    assert ".vcf" in reply


def test_photo_with_disallowed_extension_is_dropped_silently(context):
    update = _update(document=_document("evil.exe", "image/png"))
    asyncio.run(roboterri.handle_photo(update, context))
    update.message.reply_text.assert_not_awaited()
