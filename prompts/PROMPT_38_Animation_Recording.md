# Prompt #38: Animation Recording and Export

**Phase:** 2.75D (Week 2)  
**Effort:** ~7 hours  
**Priority:** Medium (publication/presentation)  
**Status:** Ready to implement  
**Dependencies:** Prompt #37 viewport rendering in place; Phase 2.75A-C complete

---

## Objective

Record simulations to GIF, MP4, or PNG sequences with configurable framerate (10-60 fps), resolution scaling (0.5x-2.0x), progress feedback, and optional auto-play. Use FFmpeg for MP4 when available; fall back gracefully to GIF/PNG when absent.

---

## Data Structures

```cpp
enum RecordingFormat { RECORDING_GIF, RECORDING_MP4, RECORDING_PNG_SEQUENCE };
enum RecordingState { RECORDING_IDLE, RECORDING_ACTIVE, RECORDING_PROCESSING, RECORDING_ERROR };

struct AnimationRecorder {
    RecordingState state;
    RecordingFormat format;
    float framerate;          // 10-60
    int target_frames;        // duration * framerate
    int frame_count;
    float resolution_scale;   // 0.5x - 2.0x
    std::vector<SDL_Surface*> frames;
    int capture_width;
    int capture_height;
    char output_path[260];
    bool auto_play;
    float progress;           // 0-1
    char status_message[128];
};

struct AppState {
    AnimationRecorder recorder;
    bool show_recording_panel;
};
```

---

## Capture Flow

1) User opens Tools -> Animation Recording panel.  
2) User selects format (GIF/MP4/PNG), framerate, duration, resolution scale, auto-play.  
3) Start recording: initialize recorder, clear previous frames, pre-compute target frame count and output path.  
4) During render loop, capture frames at the target interval using `SDL_RenderReadPixels` and optional scaling to even dimensions.  
5) When target frames reached or Stop clicked, transition to PROCESSING and encode:
   - GIF: jo_gif or stb_image_write
   - MP4: write PNG temp sequence then call FFmpeg
   - PNG: direct numbered sequence  
6) Update status/progress; on success optionally auto-open the result; on failure set error state and message.

---

## Capture Helper (Pseudo-code)

```cpp
static SDL_Surface* capture_frame(RenderContext* render, const AnimationRecorder* rec) {
    int w, h;
    SDL_GetRendererOutputSize(render->renderer, &w, &h);

    int out_w = (int)(w * rec->resolution_scale);
    int out_h = (int)(h * rec->resolution_scale);
    out_w = (out_w / 2) * 2;
    out_h = (out_h / 2) * 2;

    SDL_Surface* surface = SDL_CreateRGBSurface(0, out_w, out_h, 24, 0x00FF0000, 0x0000FF00, 0x000000FF, 0);
    if (!surface) return nullptr;

    SDL_Rect src = {0, 0, w, h};
    if (SDL_RenderReadPixels(render->renderer, &src, SDL_PIXELFORMAT_RGB24, surface->pixels, surface->pitch) != 0) {
        SDL_FreeSurface(surface);
        return nullptr;
    }

    if (out_w != w || out_h != h) {
        SDL_Surface* scaled = SDL_CreateRGBSurface(0, out_w, out_h, 24, 0x00FF0000, 0x0000FF00, 0x000000FF, 0);
        if (scaled) {
            SDL_BlitScaled(surface, nullptr, scaled, nullptr);
            SDL_FreeSurface(surface);
            return scaled;
        }
    }
    return surface;
}
```

---

## Start / Capture / Stop (Pseudo-code)

```cpp
static void start_recording(AppState* app, RecordingFormat format, float fps, int duration_sec, float scale) {
    AnimationRecorder* rec = &app->recorder;
    rec->state = RECORDING_ACTIVE;
    rec->format = format;
    rec->framerate = fps;
    rec->target_frames = (int)(fps * duration_sec);
    rec->frame_count = 0;
    rec->resolution_scale = scale;
    rec->frames.clear();
    rec->frames.reserve(rec->target_frames);

    time_t now = time(nullptr);
    struct tm* tm_info = localtime(&now);
    const char* ext = (format == RECORDING_GIF) ? "gif" : (format == RECORDING_MP4) ? "mp4" : "png";
    snprintf(rec->output_path, sizeof(rec->output_path),
             "recordings/emwave_%04d%02d%02d_%02d%02d%02d.%s",
             tm_info->tm_year + 1900, tm_info->tm_mon + 1, tm_info->tm_mday,
             tm_info->tm_hour, tm_info->tm_min, tm_info->tm_sec, ext);
    snprintf(rec->status_message, sizeof(rec->status_message), "Recording: 0/%d frames", rec->target_frames);
}

static void capture_frame_if_recording(AppState* app, RenderContext* render) {
    AnimationRecorder* rec = &app->recorder;
    if (rec->state != RECORDING_ACTIVE) return;

    SDL_Surface* frame = capture_frame(render, rec);
    if (!frame) {
        rec->state = RECORDING_ERROR;
        snprintf(rec->status_message, sizeof(rec->status_message), "Failed to capture frame");
        return;
    }
    rec->frames.push_back(frame);
    rec->frame_count++;
    rec->progress = (float)rec->frame_count / (float)rec->target_frames;
    snprintf(rec->status_message, sizeof(rec->status_message), "Recording: %d/%d (%.0f%%)",
             rec->frame_count, rec->target_frames, rec->progress * 100.0f);

    if (rec->frame_count >= rec->target_frames) stop_recording(app);
}
```

---

## Export Helpers (Outline)

- **GIF:** Include `jo_gif.h` (single header). Write all frames with delay = 1000/framerate ms.  
- **MP4:** Write PNG temp sequence to `recordings/temp/frame_%04d.png`, then call  
  `ffmpeg -y -framerate <fps> -i recordings/temp/frame_%04d.png -c:v libx264 -pix_fmt yuv420p "<output_path>"`.  
  Remove temp files after encode; set state to ERROR if return code != 0.  
- **PNG Sequence:** Write numbered PNG files to `<output_path>_frames/`.  
- Free all `SDL_Surface*` after encoding.

---

## Recording Panel (ImGui)

- Sliders: Duration (1-60 s), Framerate (10-60 fps), Resolution scale (0.5x-2.0x)  
- Combo: Format (GIF, MP4 [requires FFmpeg], PNG Sequence)  
- Checkbox: Auto-play after export  
- Buttons: Start Recording, Stop Recording  
- Progress bar + status text for ACTIVE/PROCESSING states  
- Error message styling when state == RECORDING_ERROR  

---

## Testing Checklist (12 cases)

- [ ] Panel opens from Tools menu and closes without side effects  
- [ ] Duration slider updates target frame count correctly  
- [ ] Framerate slider respected (verify capture timing)  
- [ ] Resolution scale applies to output dimensions (even numbers)  
- [ ] GIF export produces playable file for 10s @ 30 fps  
- [ ] MP4 export succeeds when FFmpeg is on PATH  
- [ ] MP4 export falls back to GIF/PNG with clear message when FFmpeg missing  
- [ ] PNG sequence writes all frames with correct naming  
- [ ] Progress bar reflects capture progress and processing phase  
- [ ] Stop button halts recording early and still saves what was captured  
- [ ] Auto-play opens output after encode (or no-op if disabled)  
- [ ] No memory leaks or crashes after freeing frames

---

## Success Criteria

1. Record/Stop controls work reliably for GIF, MP4, and PNG sequence.  
2. Framerate, duration, and resolution settings are applied to outputs.  
3. FFmpeg absence does not block recording; fallback path succeeds.  
4. Progress and status messaging are accurate and user-friendly.  
5. Encoded files open in standard players/editors without post-fixes.  
6. Minimal performance impact while recording; no persistent memory growth.

*Estimated effort: 7 hours*  
