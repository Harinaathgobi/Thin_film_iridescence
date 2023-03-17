#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_out = dot(half_vector, dir_out);
    if (n_dot_out <= 0 || n_dot_h <= 0 || h_dot_out <= 0) {
        return make_zero_spectrum();
    }

    // They use a hardcoded IOR 1.5 -> F0 = 0.04
    Real F = schlick_fresnel(Real(0.04), h_dot_out);
    // Generalized Trowbridge-Reitz distribution
    Real D = GTR1(n_dot_h, alpha);
    // SmithG with fixed alpha
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), Real(0.25), Real(0.25));
    Real G_out = smith_masking_gtr2(to_local(frame, dir_out), Real(0.25), Real(0.25));
    Real G = G_in * G_out;

    return make_const_spectrum(F * D * G / (4 * dot(dir_in, vertex.shading_frame.n)));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_out = dot(half_vector, dir_out);
    if (n_dot_out <= 0 || n_dot_h <= 0 || h_dot_out <= 0) {
        return 0;
    }

    // We only importance sample D
    Real D = GTR1(n_dot_h, alpha);

    return D * n_dot_h / (4 * h_dot_out);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    // Since the clearcoat BRDF does not have a clear geometry meaning,
    // it is hard to sample the visible normals. We will just importance sample D
    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    // Appendix B.2 Burley's note
    Real alpha2 = alpha * alpha;
    // Equation 5
    Real cos_h_elevation =
        sqrt(max(Real(0), (1 - pow(alpha2, 1 - rnd_param_uv[0])) / (1 - alpha2)));
    Real sin_h_elevation = sqrt(max(1 - cos_h_elevation * cos_h_elevation, Real(0)));
    Real h_azimuth = 2 * c_PI * rnd_param_uv[1];
    Vector3 local_micro_normal{
        sin_h_elevation * cos(h_azimuth),
        sin_h_elevation * sin(h_azimuth),
        cos_h_elevation
    };
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);

    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected, Real(0) /* eta */, sqrt(alpha) /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
