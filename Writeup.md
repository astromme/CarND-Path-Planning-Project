Notes
- Started by having the car drive straight
- Then tried curve
- Then tried following waypoints
- Then tried interpolating waypoints



Pseudocode:

- State machine [straight, pass_right, pass_left]
- Go straight at constant speed
-- unless car is in the way.
-- if car is in the way, follow it at safe distance
--- when safe to change lanes, do that
