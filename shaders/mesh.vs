#version 330 core
layout (location = 0) in vec3 pos;

uniform mat4 world2screen;

// after transformation
// x: [-0.9, 0.9]
// y:[-0.9, 0.9]
// z/h: [0, 1]

void main()
{
    gl_Position = world2screen*vec4(pos, 1.0);
}