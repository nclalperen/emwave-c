# Prompt #34: Grid Overlay Enhancements

**Phase:** 2.75C (Polish - Week 3)
**Effort:** ~3 hours
**Priority:** LOW - Visual polish
**Status:** Ready for implementation
**Dependencies:** Phase 2.75A-B complete

---

## Objective

Enhance grid overlay with adaptive density, coordinate labels, rulers, and origin marker for professional CAD-like appearance.

---

## Features

### 1. Adaptive Grid Density

```cpp
// Adjust grid spacing based on zoom level
float grid_spacing = calculate_grid_spacing(scale, sim->lx, sim->nx);

// At low zoom (0.5×-1×): Sparse grid (every 10 cells)
// At medium zoom (2×-4×): Normal grid (every 5 cells)
// At high zoom (8×+): Dense grid (every cell)

int grid_step = (scale < 2) ? 10 : (scale < 8) ? 5 : 1;
```

### 2. Coordinate Labels

```cpp
// Draw labels at grid intersections
for (int i = 0; i < sim->nx; i += grid_step) {
    float x_meters = (float)i * (sim->lx / sim->nx);
    char label[32];
    snprintf(label, sizeof(label), "%.3f m", x_meters);

    // Draw text at grid line
    ImGui::SetCursorPos(...);
    ImGui::Text("%s", label);
}
```

### 3. Rulers (Top & Left Edges)

```cpp
// Top ruler
ImDrawList* draw_list = ImGui::GetWindowDrawList();
for (int i = 0; i < sim->nx; i += grid_step) {
    float x_px = i * scale + viewport_offset_x;
    // Draw tick mark
    draw_list->AddLine(ImVec2(x_px, 0), ImVec2(x_px, 10), color);
    // Draw label
}

// Left ruler (similar)
```

### 4. Origin Marker

```cpp
// Draw cross at (0, 0)
float origin_x = viewport_offset_x;
float origin_y = viewport_offset_y;

draw_list->AddLine(ImVec2(origin_x - 10, origin_y),
                   ImVec2(origin_x + 10, origin_y),
                   IM_COL32(255, 0, 0, 255), 2.0f);  // Red cross
draw_list->AddLine(ImVec2(origin_x, origin_y - 10),
                   ImVec2(origin_x, origin_y + 10),
                   IM_COL32(255, 0, 0, 255), 2.0f);
```

---

## Testing

- [ ] Grid density adapts to zoom level
- [ ] Labels readable at all zoom levels
- [ ] Rulers show correct measurements
- [ ] Origin marker visible
- [ ] No performance impact (<1ms)

---

## Success Criteria

1. ✅ Adaptive grid spacing
2. ✅ Coordinate labels on grid
3. ✅ Rulers on edges
4. ✅ Origin cross marker
5. ✅ Professional CAD appearance

---

*Created: 2025-11-20*
