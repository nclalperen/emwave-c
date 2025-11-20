# Prompt #39: Advanced Measurement Tools

**Phase:** 2.75D (Week 2)  
**Effort:** ~5 hours  
**Priority:** Medium (research/documentation)  
**Status:** Ready to implement  
**Dependencies:** Prompt #35 ruler baseline; Prompt #37 viewport focus for correct hit-tests

---

## Objective

Extend measurements with polygon area (shoelace), perimeter calculation, text annotations with custom colors, a history panel, and CSV export for analysis.

---

## Data Structures

```cpp
enum MeasurementType { MEASURE_DISTANCE, MEASURE_AREA, MEASURE_PERIMETER, MEASURE_ANNOTATION };

struct DistanceMeasurement { ImVec2 a, b; double distance_m; double angle_deg; };
struct AreaMeasurement { std::vector<ImVec2> vertices; bool closed; double area_m2; double perimeter_m; };
struct Annotation { ImVec2 grid_pos; char text[128]; ImVec4 color; float font_size; bool visible; };

struct MeasurementHistory {
    std::vector<DistanceMeasurement> distances;
    std::vector<AreaMeasurement> areas;
    std::vector<Annotation> annotations;
};

struct AppState {
    bool area_mode;
    AreaMeasurement current_area;
    bool annotation_mode;
    Annotation temp_annotation;
    MeasurementHistory measurements;
    bool show_measurement_history;
    bool show_distance_measurements;
    bool show_area_measurements;
    bool show_annotations;
};
```

---

## Area Tool (Shoelace)

```cpp
// Toggle with 'A'
case SDLK_a: app.area_mode = !app.area_mode; break;

// Add vertex on left-click inside active viewport
if (app.area_mode && left_click && viewport_hit) {
    ImVec2 grid = screen_to_grid(mouse, scale, offset);
    app.current_area.vertices.push_back(grid);
    if (app.current_area.vertices.size() >= 3) {
        float dist = distance_to_first_vertex(app.current_area, grid);
        if (dist < 3.0f) close_area_measurement(&app, sim);
    }
}

// Right-click closes polygon if >=3 vertices
if (app.area_mode && right_click && app.current_area.vertices.size() >= 3) {
    close_area_measurement(&app, sim);
}

static void close_area_measurement(AppState* app, const SimulationState* sim) {
    AreaMeasurement& a = app->current_area;
    if (a.vertices.size() < 3) return;
    a.closed = true;

    double sum = 0.0;
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        size_t j = (i + 1) % a.vertices.size();
        sum += a.vertices[i].x * a.vertices[j].y;
        sum -= a.vertices[j].x * a.vertices[i].y;
    }
    double area_cells = std::abs(sum) * 0.5;
    double dx = sim->lx / (double)sim->nx;
    double dy = sim->ly / (double)sim->ny;
    double cell_area_m2 = dx * dy;
    a.area_m2 = area_cells * cell_area_m2;

    double per_cells = 0.0;
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        size_t j = (i + 1) % a.vertices.size();
        double dx_cells = a.vertices[j].x - a.vertices[i].x;
        double dy_cells = a.vertices[j].y - a.vertices[i].y;
        per_cells += std::sqrt(dx_cells * dx_cells + dy_cells * dy_cells);
    }
    double avg_cell_m = 0.5 * (dx + dy);
    a.perimeter_m = per_cells * avg_cell_m;

    app->measurements.areas.push_back(a);
    a.vertices.clear();
    a.closed = false;
    app->area_mode = false;
}
```

---

## Drawing Overlays

```cpp
// Convert grid to screen for the active viewport
auto grid_to_screen = [&](ImVec2 grid) -> ImVec2 {
    ImVec2 offset = compute_viewport_offset(app, sim, scale);
    return ImVec2(app.viewport_pos.x + offset.x + grid.x * (float)scale,
                  app.viewport_pos.y + offset.y + grid.y * (float)scale);
};

// Saved areas
for (const AreaMeasurement& area : app.measurements.areas) {
    if (area.vertices.size() < 3) continue;
    std::vector<ImVec2> pts;
    for (ImVec2 v : area.vertices) pts.push_back(grid_to_screen(v));
    dl->AddConvexPolyFilled(pts.data(), (int)pts.size(), IM_COL32(80, 180, 80, 80));
    for (size_t i = 0; i < pts.size(); ++i) {
        size_t j = (i + 1) % pts.size();
        dl->AddLine(pts[i], pts[j], IM_COL32(100, 255, 100, 200), 2.0f);
    }
    // Label at centroid
    ImVec2 centroid = compute_centroid(pts);
    char buf[96];
    snprintf(buf, sizeof(buf), "%.6f m2\n%.4f m", area.area_m2, area.perimeter_m);
    draw_label(dl, centroid, buf, IM_COL32(0,0,0,200), IM_COL32(120,255,120,255));
}

// In-progress polygon (when area_mode is true) with preview line to cursor
```

---

## Text Annotations

```cpp
// Activate with 'T'
case SDLK_t:
    app.annotation_mode = true;
    app.temp_annotation = {};
    app.temp_annotation.color = ImVec4(1, 1, 0.2f, 1);
    app.temp_annotation.font_size = 14.0f;
    app.temp_annotation.grid_pos = current_grid_cursor;
    ImGui::OpenPopup("AddAnnotation");
    break;

if (ImGui::BeginPopup("AddAnnotation")) {
    ImGui::InputText("Text", app.temp_annotation.text, sizeof(app.temp_annotation.text));
    ImGui::ColorEdit3("Color", (float*)&app.temp_annotation.color);
    ImGui::SliderFloat("Size", &app.temp_annotation.font_size, 8.0f, 28.0f, "%.0f px");
    if (ImGui::Button("Add")) { app.temp_annotation.visible = true; app.measurements.annotations.push_back(app.temp_annotation); ImGui::CloseCurrentPopup(); app.annotation_mode = false; }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) { ImGui::CloseCurrentPopup(); app.annotation_mode = false; }
    ImGui::EndPopup();
}
```

Drawing:
```cpp
for (const Annotation& ann : app.measurements.annotations) {
    if (!ann.visible) continue;
    ImVec2 screen = grid_to_screen(ann.grid_pos);
    ImVec2 size = ImGui::CalcTextSize(ann.text);
    ImVec2 pad(6, 3);
    ImVec2 min = ImVec2(screen.x - pad.x, screen.y - pad.y);
    ImVec2 max = ImVec2(screen.x + size.x + pad.x, screen.y + size.y + pad.y);
    dl->AddRectFilled(min, max, IM_COL32(0,0,0,200), 3.0f);
    dl->AddRect(min, max, ImColor(ann.color), 3.0f, 0, 2.0f);
    dl->AddText(ImVec2(screen.x, screen.y), ImColor(ann.color), ann.text);
}
```

---

## Measurement History Panel

- Checkboxes: Show Distances, Show Areas, Show Annotations  
- Sections:
  - Distances: list with delete buttons; values shown with distance/angle.  
  - Areas: list with area/perimeter; delete button.  
  - Annotations: visibility toggles + delete button.  
- Buttons: Export to CSV, Clear All  

CSV export (outline):
```cpp
static void export_measurements_csv(const AppState* app, const char* path) {
    FILE* f = fopen(path, "w");
    fprintf(f, "Type,ID,Value,Unit,Data\n");
    for (size_t i = 0; i < app->measurements.distances.size(); ++i) {
        const auto& d = app->measurements.distances[i];
        fprintf(f, "Distance,%zu,%.6f,m,\"(%.1f,%.1f)-(%.1f,%.1f) @ %.1f deg\"\n",
                i+1, d.distance_m, d.a.x, d.a.y, d.b.x, d.b.y, d.angle_deg);
    }
    for (size_t i = 0; i < app->measurements.areas.size(); ++i) {
        const auto& a = app->measurements.areas[i];
        fprintf(f, "Area,%zu,%.6f,m2,\"%zu vertices | Perimeter=%.4f m\"\n",
                i+1, a.area_m2, a.vertices.size(), a.perimeter_m);
    }
    for (size_t i = 0; i < app->measurements.annotations.size(); ++i) {
        const auto& ann = app->measurements.annotations[i];
        fprintf(f, "Annotation,%zu,,,\"(%.1f,%.1f): %s\"\n",
                i+1, ann.grid_pos.x, ann.grid_pos.y, ann.text);
    }
    fclose(f);
}
```

---

## Testing Checklist (15 cases)

- [ ] 'A' toggles area mode on/off  
- [ ] Left-click adds vertices only inside active viewport  
- [ ] Right-click or closing near first vertex finalizes polygon  
- [ ] Area computed correctly on rectangle (compare expected)  
- [ ] Perimeter computed correctly on rectangle/triangle  
- [ ] Area fill, outline, and vertex markers render without artifacts  
- [ ] Live preview line shows from last vertex to cursor while drawing  
- [ ] Annotation dialog opens with 'T' and adds colored text at cursor grid position  
- [ ] Annotation visibility toggle works in history panel  
- [ ] Measurement history lists distances, areas, annotations with delete controls  
- [ ] CSV export opens in Excel/LibreOffice with correct headers/values  
- [ ] Distances from ruler tool are added to history  
- [ ] Show/hide toggles for areas/annotations respond immediately  
- [ ] Handles 50+ measurements without frame drops  
- [ ] Clearing all measurements removes overlays

---

## Success Criteria

1. Accurate area and perimeter outputs using grid spacing.  
2. Reliable annotation creation with configurable color and size.  
3. History panel allows review/delete of all measurement types.  
4. CSV export produces clean, analysis-ready data.  
5. Overlays render correctly across all viewport layouts from Prompt #37.  
6. Performance remains smooth with dense measurement sets.

*Estimated effort: 5 hours*  
