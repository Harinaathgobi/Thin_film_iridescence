#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Common variables
    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_out = dot(half_vector, dir_out);
    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);

    // Disney's metallic BRDF is standard Cook-Torrence
    // microfacet model: F * G * D / (4 * n_dot_in)
    // (We include n_dot_out in our BRDFs)

    // For F, they use a Schlick approximation
    // F0 + (1 - F0) * (1 - h_dot_out)
    // For F0 they blend between base color and RGB(1, 1, 1) using
    // specular and metallic parameters. This does not make much sense
    // in our case since we are only modelling the metallic part --
    // we will use just the base color.
    Spectrum F = schlick_fresnel(base_color, h_dot_out);

    // For D, they use an anisotropic Trowbridge-Reitz distribution (aka GGX)
    Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);

    // For G, they use the Smith masking term corresponds to GTR2
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
    Real G_out = smith_masking_gtr2(to_local(frame, dir_out), alpha_x, alpha_y);
    Real G = G_in * G_out;

    // And we're done!
    return F * D * G / (4 * dot(dir_in, vertex.shading_frame.n));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Common variables
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);

    // We use visible normal sampling, so the PDF ~ (G_in * D) / (4 * n_dot_in)
    Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
    return D * G_in / (4 * dot(dir_in, vertex.shading_frame.n));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);

    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected, Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
