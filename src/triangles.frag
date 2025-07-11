#version 450 core

in vec4 color;
in vec3 fragNormal;  
in vec3 fragPosition;  
in vec3 baseColor; 

uniform vec3 lightDir;
uniform int shadingMode;

out vec4 fColor;

void main()
{
    if (shadingMode == 2) { // PHONG
        // Normalize the interpolated normal
        vec3 normal = normalize(fragNormal);
        
        // Ambient
        float ambientStrength = 0.2;
        vec3 ambient = ambientStrength * baseColor;
        
        // Diffuse
        float diff = max(dot(normal, lightDir), 0.0);
        vec3 diffuse = diff * baseColor;
        
        // Specular (Blinn-Phong)
        float specularStrength = 0.5;
        vec3 viewDir = normalize(-fragPosition); // In view space, camera is at (0,0,0)
        vec3 halfwayDir = normalize(lightDir + viewDir);
        float spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);
        vec3 specular = specularStrength * spec * vec3(1.0);
        
        vec3 result = ambient + diffuse + specular;
        fColor = vec4(result, 1.0);
    }
    else {
        // For NO_SHADING and GOURAUD, just use the color calculated in vertex shader
        fColor = color;
    }
}
