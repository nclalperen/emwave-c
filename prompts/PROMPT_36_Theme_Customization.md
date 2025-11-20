# Prompt #36: Theme & Style Customization

**Phase:** 2.75C (Polish - Week 3 - Final)
**Effort:** ~4 hours
**Priority:** LOW - Visual polish
**Status:** Ready for implementation
**Dependencies:** Phase 2.75A-B complete

---

## Objective

Enhance UI themes with rounded corners, larger fonts (4K support), and Blender-inspired colors for professional appearance.

---

## Implementation

### 1. Enhanced Theme System

```cpp
enum ThemePreset {
    THEME_DARK_PROFESSIONAL = 0,
    THEME_BLENDER_INSPIRED = 1,
    THEME_LIGHT_MODE = 2,
    THEME_HIGH_CONTRAST = 3
};

static void apply_theme(ThemePreset theme) {
    ImGuiStyle& style = ImGui::GetStyle();

    // Rounded corners (modern look)
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.GrabRounding = 3.0f;
    style.TabRounding = 3.0f;

    // Padding and spacing
    style.WindowPadding = ImVec2(8, 8);
    style.FramePadding = ImVec2(8, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.ItemInnerSpacing = ImVec2(6, 6);

    ImVec4* colors = style.Colors;

    switch (theme) {
        case THEME_DARK_PROFESSIONAL:
            colors[ImGuiCol_WindowBg] = ImVec4(0.10f, 0.10f, 0.11f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.12f, 0.12f, 0.13f, 1.00f);
            colors[ImGuiCol_FrameBg] = ImVec4(0.18f, 0.18f, 0.20f, 1.00f);
            colors[ImGuiCol_FrameBgHovered] = ImVec4(0.25f, 0.25f, 0.28f, 1.00f);
            colors[ImGuiCol_FrameBgActive] = ImVec4(0.30f, 0.30f, 0.35f, 1.00f);
            colors[ImGuiCol_TitleBg] = ImVec4(0.08f, 0.08f, 0.09f, 1.00f);
            colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.12f, 0.15f, 1.00f);
            colors[ImGuiCol_Button] = ImVec4(0.20f, 0.20f, 0.25f, 1.00f);
            colors[ImGuiCol_ButtonHovered] = ImVec4(0.28f, 0.28f, 0.35f, 1.00f);
            colors[ImGuiCol_ButtonActive] = ImVec4(0.35f, 0.35f, 0.45f, 1.00f);
            colors[ImGuiCol_Header] = ImVec4(0.22f, 0.22f, 0.28f, 1.00f);
            colors[ImGuiCol_HeaderHovered] = ImVec4(0.28f, 0.28f, 0.35f, 1.00f);
            colors[ImGuiCol_HeaderActive] = ImVec4(0.32f, 0.32f, 0.40f, 1.00f);
            break;

        case THEME_BLENDER_INSPIRED:
            // Blender 3.x colors (dark gray + blue accent)
            colors[ImGuiCol_WindowBg] = ImVec4(0.16f, 0.16f, 0.16f, 1.00f);
            colors[ImGuiCol_ChildBg] = ImVec4(0.18f, 0.18f, 0.18f, 1.00f);
            colors[ImGuiCol_TitleBgActive] = ImVec4(0.26f, 0.46f, 0.76f, 1.00f);
            colors[ImGuiCol_Button] = ImVec4(0.28f, 0.28f, 0.28f, 1.00f);
            colors[ImGuiCol_ButtonHovered] = ImVec4(0.38f, 0.58f, 0.88f, 1.00f);
            colors[ImGuiCol_ButtonActive] = ImVec4(0.26f, 0.46f, 0.76f, 1.00f);
            // ... more Blender colors
            break;

        case THEME_LIGHT_MODE:
            ImGui::StyleColorsLight();
            // Adjust for better readability
            colors[ImGuiCol_WindowBg] = ImVec4(0.95f, 0.95f, 0.95f, 1.00f);
            break;

        case THEME_HIGH_CONTRAST:
            // Black background, white text, high saturation
            colors[ImGuiCol_WindowBg] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
            colors[ImGuiCol_Text] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
            colors[ImGuiCol_Button] = ImVec4(0.2f, 0.2f, 0.8f, 1.0f);
            break;
    }
}
```

---

### 2. Font Scaling for 4K Displays

```cpp
// In initialization
ImGuiIO& io = ImGui::GetIO();

// Detect DPI scaling
float dpi_scale = 1.0f;
int display_index = SDL_GetWindowDisplayIndex(window);
float ddpi, hdpi, vdpi;
if (SDL_GetDisplayDPI(display_index, &ddpi, &hdpi, &vdpi) == 0) {
    dpi_scale = ddpi / 96.0f;  // 96 DPI = 1.0×
}

// Load larger font for high DPI
if (dpi_scale > 1.5f) {
    io.Fonts->AddFontFromFileTTF("assets/fonts/DejaVuSans.ttf", 18.0f * dpi_scale);
} else {
    io.Fonts->AddFontDefault();
}

io.FontGlobalScale = 1.0f;  // Fine-tune if needed
```

---

### 3. Theme Menu

```cpp
// In menu bar
if (ImGui::BeginMenu("View")) {
    if (ImGui::BeginMenu("Theme")) {
        if (ImGui::MenuItem("Dark Professional")) {
            apply_theme(THEME_DARK_PROFESSIONAL);
        }
        if (ImGui::MenuItem("Blender Inspired")) {
            apply_theme(THEME_BLENDER_INSPIRED);
        }
        if (ImGui::MenuItem("Light Mode")) {
            apply_theme(THEME_LIGHT_MODE);
        }
        if (ImGui::MenuItem("High Contrast")) {
            apply_theme(THEME_HIGH_CONTRAST);
        }
        ImGui::EndMenu();
    }
    ImGui::EndMenu();
}
```

---

### 4. Color Accent Customization

```cpp
// AppState addition
struct AppState {
    ImVec4 accent_color;  // User-customizable accent
};

// Apply accent to key elements
colors[ImGuiCol_ButtonHovered] = app.accent_color;
colors[ImGuiCol_HeaderHovered] = app.accent_color;
colors[ImGuiCol_TabActive] = app.accent_color;

// Color picker in settings
ImGui::ColorEdit3("Accent Color", (float*)&app.accent_color);
```

---

## Testing Checklist

- [ ] Dark Professional theme - clean, modern
- [ ] Blender Inspired theme - matches Blender 3.x
- [ ] Light Mode theme - readable on bright displays
- [ ] High Contrast theme - accessibility
- [ ] Rounded corners visible on panels
- [ ] Font scaling works on 4K displays
- [ ] Theme switcher in View menu
- [ ] Accent color customizable
- [ ] All themes readable and professional

---

## Success Criteria

1. ✅ 4 theme presets implemented
2. ✅ Rounded corners (modern style)
3. ✅ 4K DPI scaling support
4. ✅ Blender-inspired theme
5. ✅ Accent color customization
6. ✅ Professional, polished appearance

---

**🎉 Phase 2.75C COMPLETE!**
**🎉 UI MODERNIZATION COMPLETE!**

**Total: Prompts #28-36 (9 prompts, ~44 hours)**

---

*Created: 2025-11-20*
