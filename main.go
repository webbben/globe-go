package globego

/*
#cgo LDFLAGS: -L./bin -lglobe
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Settings struct for configuring Globe
typedef struct {
    size_t refresh_rate;
    float globe_rotation_speed;
    float cam_rotation_speed;
    float cam_zoom;
    float focus_speed;
    bool night;
    float coord_x;
	float coord_y;
} SettingsFFI;

// Declare the Rust functions here
void ext_screensaver(SettingsFFI settings);
void ext_interactive(SettingsFFI settings);
*/
import "C"

func Screensaver(settings C.SettingsFFI) {
	C.ext_screensaver(settings)
}

func Interactive(settings C.SettingsFFI) {
	C.ext_interactive(settings)
}

func GetDefaultSettings() C.SettingsFFI {
	return C.SettingsFFI{
		refresh_rate:         60,
		globe_rotation_speed: 0.01,
		cam_rotation_speed:   0,
		cam_zoom:             0.1,
		focus_speed:          1,
		night:                false,
		coord_x:              35,
		coord_y:              135,
	}
}
