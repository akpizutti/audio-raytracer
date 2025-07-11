
#version 400 core

layout( location = 0 ) in vec4 vPosition;
layout( location = 1 ) in vec3 vNormal;



//out vec4 position;
//out vec4 normalVertex;
out vec4 color;
out vec3 fragNormal;   
out vec3 fragPosition;  
out vec3 baseColor; 

uniform mat4 projection;
uniform mat4 model;
uniform mat4 view;
uniform vec3 in_color;
uniform vec3 lightDir;
uniform int shadingMode; 


void main()
{
    vec4 worldPos = model * vPosition;
    gl_Position = projection * view * worldPos;
    
    // Transform normal to view space
    mat3 normalMatrix = transpose(inverse(mat3(view * model)));
    vec3 normal = normalize(normalMatrix * vNormal);
    
    // Pass variables needed for Phong shading to the fragment shader
    fragNormal = normal;
    fragPosition = vec3(view * worldPos);
    baseColor = in_color;
    
    // Shading modes:
    // 0 = NO_SHADING - just use the base color
    // 1 = GOURAUD - compute lighting in vertex shader
    // 2 = PHONG - compute lighting in fragment shader
    
    if (shadingMode == 0) { // NO_SHADING
        color = vec4(in_color, 1.0);
    }
    else if (shadingMode == 1) { // GOURAUD
        // Calculate lighting in the vertex shader
        float ambientStrength = 0.2;
        vec3 ambient = ambientStrength * in_color;
        
        // Diffuse
        float diff = max(dot(normal, lightDir), 0.0);
        vec3 diffuse = diff * in_color;
        
        // Specular (simple Blinn-Phong)
        float specularStrength = 0.5;
        vec3 viewDir = vec3(0.0, 0.0, 1.0); // In view space, camera is at (0,0,0) looking along -z
        vec3 halfwayDir = normalize(lightDir + viewDir);
        float spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);
        vec3 specular = specularStrength * spec * vec3(1.0);
        
        vec3 result = ambient + diffuse + specular;
        color = vec4(result, 1.0);
    }
    else { // PHONG
        // For Phong shading, just pass data to the fragment shader
        // The actual lighting calculation will be done in the fragment shader
        color = vec4(0.0); // Will be computed in fragment shader
    }
}