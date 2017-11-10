## Path Planning

At a high level, paths are created by interpolating between the waypoints using
a spline, adjusting the lane as necessary to maximize the distance between the
car and cars ahead of it.

In depth:

1. I keep track of a speed, and increase it by an acceleration factor (0.125) at each timestep until the speed matches the target speed. If there is a car ahead, I reduce the speed until that car is further than 30 meters. This basic behavior causes some oscillation because the car overshoots the ideal speed because distance doesn't change instantaneously with speed.

2. This speed is used to plot the next point by computing the delta between the current point and the next point. A local reference frame is used where the forward direction is x. A spline is used to get y from x.

3. The spline is computed by starting with two waypoints that are tangential to the car's current yaw, and then adding on a few more waypoints in the future. These future waypoints don't have to be in the same lane, and if they aren't, the spline naturally interpolates gradually.

4. The model has a cooldown of 50 timesteps between lane changes to prevent changing too fast.

5. The ideal lane is computed by observing which lane has the most open road without another car. The car is greedy at switching lanes as long as there is not a car within 30 meters ahead in the new lane.

6. To handle the edge case of the end of the map, I add a few extra waypoints with increasing s but wrapped around x, y. This lets the car see into the future even across the boundary where the car is at high s and the waypoints in the future are just past zero s.
