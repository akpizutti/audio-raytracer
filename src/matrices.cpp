#include "matrices.h"

// Esta fun��o Matrix() auxilia na cria��o de matrizes usando a biblioteca GLM.
// Note que em OpenGL (e GLM) as matrizes s�o definidas como "column-major",
// onde os elementos da matriz s�o armazenadas percorrendo as COLUNAS da mesma.
// Por exemplo, seja
//
//       [a b c]
//   M = [d e f]
//       [g h i]
//
// uma matriz 3x3. Em mem�ria, na representa��o "column-major" de OpenGL, essa
// matriz � representada pelo seguinte array:
//
//   M[] = {  a,d,g,    b,e,h,    c,f,i  };
//              ^         ^         ^
//              |         |         |
//           coluna 1  coluna 2  coluna 3
//
// Para conseguirmos definir matrizes atrav�s de suas LINHAS, a fun��o Matrix()
// computa a transposta usando os elementos passados por par�metros.
glm::mat4 Matrix(
    float m00, float m01, float m02, float m03, // LINHA 1
    float m10, float m11, float m12, float m13, // LINHA 2
    float m20, float m21, float m22, float m23, // LINHA 3
    float m30, float m31, float m32, float m33  // LINHA 4
)
{
    return glm::mat4(
        m00, m10, m20, m30, // COLUNA 1
        m01, m11, m21, m31, // COLUNA 2
        m02, m12, m22, m32, // COLUNA 3
        m03, m13, m23, m33  // COLUNA 4
    );
}

// Matriz identidade.
glm::mat4 Matrix_Identity()
{
    return Matrix(
        1.0f, 0.0f, 0.0f, 0.0f, // LINHA 1
        0.0f, 1.0f, 0.0f, 0.0f, // LINHA 2
        0.0f, 0.0f, 1.0f, 0.0f, // LINHA 3
        0.0f, 0.0f, 0.0f, 1.0f   // LINHA 4
    );
}

// Matriz de transla��o T. Seja p=[px,py,pz,pw] um ponto e t=[tx,ty,tz,0] um
// vetor em coordenadas homog�neas, definidos em um sistema de coordenadas
// Cartesiano. Ent�o, a matriz T � definida pela seguinte igualdade:
//
//     T*p = p+t.
//
glm::mat4 Matrix_Translate(float tx, float ty, float tz)
{
    return Matrix(
        1.0f, 0.0f, 0.0f, tx,
        0.0f, 1.0f, 0.0f, ty,
        0.0f, 0.0f, 1.0f, tz,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Matriz S de "escalamento de um ponto" em rela��o � origem do sistema de
// coordenadas. Seja p=[px,py,pz,pw] um ponto em coordenadas homog�neas.
// Ent�o, a matriz S � definida pela seguinte igualdade:
//
//     S*p = [sx*px, sy*py, sz*pz, pw].
//
glm::mat4 Matrix_Scale(float sx, float sy, float sz)
{
    return Matrix(
        sx, 0.0f, 0.0f, 0.0f,
        0.0f, sy, 0.0f, 0.0f,
        0.0f, 0.0f, sz, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Matriz R de "rota��o de um ponto" em rela��o � origem do sistema de
// coordenadas e em torno do eixo X (primeiro vetor da base do sistema de
// coordenadas). Seja p=[px,py,pz,pw] um ponto em coordenadas homog�neas.
// Ent�o, a matriz R � definida pela seguinte igualdade:
//
//   R*p = [ px, c*py-s*pz, s*py+c*pz, pw ];
//
// onde 'c' e 's' s�o o cosseno e o seno do �ngulo de rota��o, respectivamente.
glm::mat4 Matrix_Rotate_X(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return Matrix(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, c, -s, 0.0f,
        0.0f, s, c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Matriz R de "rota��o de um ponto" em rela��o � origem do sistema de
// coordenadas e em torno do eixo Y (segundo vetor da base do sistema de
// coordenadas). Seja p=[px,py,pz,pw] um ponto em coordenadas homog�neas.
// Ent�o, a matriz R � definida pela seguinte igualdade:
//
//   R*p = [ c*px+s*pz, py, -s*px+c*pz, pw ];
//
// onde 'c' e 's' s�o o cosseno e o seno do �ngulo de rota��o, respectivamente.
glm::mat4 Matrix_Rotate_Y(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return Matrix(
        c, 0.0f, s, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        -s, 0.0f, c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Matriz R de "rota��o de um ponto" em rela��o � origem do sistema de
// coordenadas e em torno do eixo Z (terceiro vetor da base do sistema de
// coordenadas). Seja p=[px,py,pz,pw] um ponto em coordenadas homog�neas.
// Ent�o, a matriz R � definida pela seguinte igualdade:
//
//   R*p = [ c*px-s*py, s*px+c*py, pz, pw ];
//
// onde 'c' e 's' s�o o cosseno e o seno do �ngulo de rota��o, respectivamente.
glm::mat4 Matrix_Rotate_Z(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return Matrix(
        c, -s, 0.0f, 0.0f,
        s, c, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Fun��o que calcula a norma Euclidiana de um vetor cujos coeficientes s�o
// definidos em uma base ortonormal qualquer.
float norm(glm::vec4 v)
{
    float vx = v.x;
    float vy = v.y;
    float vz = v.z;

    return sqrt(vx * vx + vy * vy + vz * vz);
}

float norm3(glm::vec3 v)
{
    float vx = v.x;
    float vy = v.y;
    float vz = v.z;

    return sqrt(vx * vx + vy * vy + vz * vz);
}

// Matriz R de "rota��o de um ponto" em rela��o � origem do sistema de
// coordenadas e em torno do eixo definido pelo vetor 'axis'. Esta matriz pode
// ser definida pela f�rmula de Rodrigues. Lembre-se que o vetor que define o
// eixo de rota��o deve ser normalizado!
glm::mat4 Matrix_Rotate(float angle, glm::vec4 axis)
{
    float c = cos(angle);
    float s = sin(angle);

    glm::vec4 v = axis / norm(axis);

    float vx = v.x;
    float vy = v.y;
    float vz = v.z;

    return Matrix(
        vx * vx * (1.0f - c) + c, vx * vy * (1.0f - c) - vz * s, vx * vz * (1 - c) + vy * s, 0.0f,
        vx * vy * (1.0f - c) + vz * s, vy * vy * (1.0f - c) + c, vy * vz * (1 - c) - vx * s, 0.0f,
        vx * vz * (1 - c) - vy * s, vy * vz * (1 - c) + vx * s, vz * vz * (1.0f - c) + c, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Produto vetorial entre dois vetores u e v definidos em um sistema de
// coordenadas ortonormal.
glm::vec4 crossproduct(glm::vec4 u, glm::vec4 v)
{
    float u1 = u.x;
    float u2 = u.y;
    float u3 = u.z;
    float v1 = v.x;
    float v2 = v.y;
    float v3 = v.z;

    return glm::vec4(
        u2 * v3 - u3 * v2, // Primeiro coeficiente
        u3 * v1 - u1 * v3, // Segundo coeficiente
        u1 * v2 - u2 * v1, // Terceiro coeficiente
        0.0f // w = 0 para vetores.
    );
}

// Produto escalar entre dois vetores u e v definidos em um sistema de
// coordenadas ortonormal.
float dotproduct(glm::vec4 u, glm::vec4 v)
{
    float u1 = u.x;
    float u2 = u.y;
    float u3 = u.z;
    float u4 = u.w;
    float v1 = v.x;
    float v2 = v.y;
    float v3 = v.z;
    float v4 = v.w;

    if (u4 != 0.0f || v4 != 0.0f)
    {
        fprintf(stderr, "ERROR: Produto escalar n�o definido para pontos.\n");
        std::exit(EXIT_FAILURE);
    }

    return u1 * v1 + u2 * v2 + u3 * v3;
}

// Matriz de mudan�a de coordenadas para o sistema de coordenadas da C�mera.
glm::mat4 Matrix_Camera_View(glm::vec4 position_c, glm::vec4 view_vector, glm::vec4 up_vector)
{
    glm::vec4 w = -view_vector;
    glm::vec4 u = crossproduct(up_vector, w);

    // Normalizamos os vetores u e w
    w = w / norm(w);
    u = u / norm(u);

    glm::vec4 v = crossproduct(w, u);

    glm::vec4 origin_o = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);

    float ux = u.x;
    float uy = u.y;
    float uz = u.z;
    float vx = v.x;
    float vy = v.y;
    float vz = v.z;
    float wx = w.x;
    float wy = w.y;
    float wz = w.z;

    return Matrix(
        ux, uy, uz, -dotproduct(u, position_c - origin_o),
        vx, vy, vz, -dotproduct(v, position_c - origin_o),
        wx, wy, wz, -dotproduct(w, position_c - origin_o),
        0.0f, 0.0f, 0.0f, 1.0f
    );
}

// Matriz de proje��o paralela ortogr�fica
glm::mat4 Matrix_Orthographic(float l, float r, float b, float t, float n, float f)
{
    glm::mat4 M = Matrix(
        2.0f / (r - l), 0.0f, 0.0f, -(r + l) / (r - l),
        0.0f, 2.0f / (t - b), 0.0f, -(t + b) / (t - b),
        0.0f, 0.0f, 2.0f / (f - n), -(f + n) / (f - n),
        0.0f, 0.0f, 0.0f, 1.0f
    );

    return M;
}

// Matriz de proje��o perspectiva
glm::mat4 Matrix_Perspective(float field_of_view, float aspect, float n, float f)
{
    float t = fabs(n) * tanf(field_of_view / 2.0f);
    float b = -t;
    float r = t * aspect;
    float l = -r;

    glm::mat4 P = Matrix(
        n, 0.0f, 0.0f, 0.0f,
        0.0f, n, 0.0f, 0.0f,
        0.0f, 0.0f, n + f, -f * n,
        0.0f, 0.0f, 1.0f, 0.0f
    );

    // A matriz M � a mesma computada acima em Matrix_Orthographic().
    glm::mat4 M = Matrix_Orthographic(l, r, b, t, n, f);

    // Note que as matrizes M*P e -M*P fazem exatamente a mesma proje��o
    // perspectiva, j� que o sinal de negativo n�o ir� afetar o resultado
    // devido � divis�o por w. Por exemplo, seja q = [qx,qy,qz,1] um ponto:
    //
    //      M*P*q = [ qx', qy', qz', w ]
    //   =(div w)=> [ qx'/w, qy'/w, qz'/w, 1 ]   Eq. (*)
    //
    // agora com o sinal de negativo:
    //
    //     -M*P*q = [ -qx', -qy', -qz', -w ]
    //   =(div w)=> [ -qx'/-w, -qy'/-w, -qz'/-w, -w/-w ]
    //            = [ qx'/w, qy'/w, qz'/w, 1 ]   Eq. (**)
    //
    // Note que o ponto final, ap�s divis�o por w, � igual: Eq. (*) == Eq. (**).
    //
    // Ent�o, por que utilizamos -M*P ao inv�s de M*P? Pois a especifica��o de
    // OpenGL define que os pontos fora do cubo unit�rio NDC dever�o ser
    // descartados j� que n�o ir�o aparecer na tela. O teste que define se um ponto
    // q est� dentro do cubo unit�rio NDC pode ser expresso como:
    //
    //      -1 <= qx'/w <= 1   &&  -1 <= qy'/w <= 1   &&  -1 <= qz'/w <= 1
    //
    // ou, de maneira equivalente SE w > 0, a placa de v�deo faz o seguinte teste
    // ANTES da divis�o por w:
    //
    //      -w <= qx' <= w   &&  -w <= qy' <= w   &&  -w <= qz' <= w
    //
    // Note que o teste acima economiza uma divis�o por w caso o ponto seja
    // descartado (quando esteja fora de NDC), entretanto, este �ltimo teste s�
    // � equivalente ao primeiro teste SE E SOMENTE SE w > 0 (isto �, se w for
    // positivo). Como este �ltimo teste � o que a placa de v�deo (GPU) ir� fazer,
    // precisamos utilizar a matriz -M*P para proje��o perspectiva, de forma que
    // w seja positivo.
    //
    return -M * P;
}

// Fun��o que imprime uma matriz M no terminal
void PrintMatrix(glm::mat4 M)
{
    printf("\n");
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ]\n", M[0][0], M[1][0], M[2][0], M[3][0]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ]\n", M[0][1], M[1][1], M[2][1], M[3][1]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ]\n", M[0][2], M[1][2], M[2][2], M[3][2]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ]\n", M[0][3], M[1][3], M[2][3], M[3][3]);
}

// Fun��o que imprime um vetor v no terminal
void PrintVector(glm::vec4 v)
{
    printf("\n");
    printf("[ %+0.2f ]\n", v[0]);
    printf("[ %+0.2f ]\n", v[1]);
    printf("[ %+0.2f ]\n", v[2]);
    printf("[ %+0.2f ]\n", v[3]);
}

// Fun��o que imprime o produto de uma matriz por um vetor no terminal
void PrintMatrixVectorProduct(glm::mat4 M, glm::vec4 v)
{
    auto r = M * v;
    printf("\n");
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]\n", M[0][0], M[1][0], M[2][0], M[3][0], v[0], r[0]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ] = [ %+0.2f ]\n", M[0][1], M[1][1], M[2][1], M[3][1], v[1], r[1]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]\n", M[0][2], M[1][2], M[2][2], M[3][2], v[2], r[2]);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]\n", M[0][3], M[1][3], M[2][3], M[3][3], v[3], r[3]);
}

// Fun��o que imprime o produto de uma matriz por um vetor, junto com divis�o
// por w, no terminal.
void PrintMatrixVectorProductDivW(glm::mat4 M, glm::vec4 v)
{
    auto r = M * v;
    auto w = r[3];
    printf("\n");
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]            [ %+0.2f ]\n", M[0][0], M[1][0], M[2][0], M[3][0], v[0], r[0], r[0] / w);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ] = [ %+0.2f ] =(div w)=> [ %+0.2f ]\n", M[0][1], M[1][1], M[2][1], M[3][1], v[1], r[1], r[1] / w);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]            [ %+0.2f ]\n", M[0][2], M[1][2], M[2][2], M[3][2], v[2], r[2], r[2] / w);
    printf("[ %+0.2f  %+0.2f  %+0.2f  %+0.2f ][ %+0.2f ]   [ %+0.2f ]            [ %+0.2f ]\n", M[0][3], M[1][3], M[2][3], M[3][3], v[3], r[3], r[3] / w);
}