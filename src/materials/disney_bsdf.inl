#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    
    // Fetch the texture values for later use
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission =
        eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = Real(0.25) * eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);
    Real alpha_c = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    if (reflect) {
        // Diffuse component
        Real n_dot_in = dot(frame.n, dir_in);
        Real n_dot_out = dot(frame.n, dir_out);
        Vector3 half_vector = normalize(dir_in + dir_out);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real n_dot_h = dot(half_vector, vertex.shading_frame.n);
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);

        Spectrum diffuse_contrib = make_zero_spectrum();
        Spectrum metallic_contrib = make_zero_spectrum();
        Real clearcoat_contrib = Real(0);
        Spectrum sheen_contrib = make_zero_spectrum();

        // For diffuse, metallic, sheen, and clearcoat, the light bounces
        // only at the top of the surface.
        if (dot(vertex.geometric_normal, dir_in) >= 0 &&
                dot(vertex.geometric_normal, dir_out) >= 0 &&
                n_dot_out > 0) {
            // Diffuse component

            // The base diffuse model
            Real Fd90 = Real(0.5) + 2 * roughness * h_dot_out * h_dot_out;
            Real schlick_n_dot_out = pow(1 - n_dot_out, Real(5));
            Real schlick_n_dot_in  = pow(1 - n_dot_in, Real(5));
            Real schlick_h_dot_out = pow(1 - h_dot_out, Real(5));
            Real base_diffuse = (1 + (Fd90 - 1) * schlick_n_dot_out) *
                                (1 + (Fd90 - 1) * schlick_n_dot_in);

            // The subsurface model
            // Disney's hack to increase the response at grazing angle
            Real Fss90 = h_dot_out * h_dot_out * roughness;
            Real Fss = (1 + (Fss90 - 1) * schlick_n_dot_out) *
                       (1 + (Fss90 - 1) * schlick_n_dot_in);
            // Lommel-Seeliger law (modified/rescaled)
            Real ss = Real(1.25) * (Fss * (1 / (n_dot_out + n_dot_in) - Real(0.5)) + Real(0.5));

            diffuse_contrib =
                (1 - specular_transmission) * (1 - metallic) * base_color *
                ((base_diffuse * (1 - subsurface) + ss * subsurface) / c_PI) * n_dot_out;

            // Sheen component
            Spectrum Ctint = 
                luminance(base_color) > 0 ? base_color / luminance(base_color) :
                fromRGB(Vector3{1, 1, 1});
            Spectrum Csheen = (1 - sheen_tint) * fromRGB(Vector3{1, 1, 1}) + sheen_tint * Ctint;
            sheen_contrib =
                (1 - metallic) * sheen * Csheen * schlick_h_dot_out * n_dot_out;

            // Metallic component
            if (n_dot_in > 0 && h_dot_out > 0 && n_dot_h > 0) {
                Real eta = bsdf.eta; // we're always going inside
                Real spec_f0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
                Spectrum spec_color =
                    ((1 - specular_tint) * fromRGB(Vector3{1, 1, 1}) + specular_tint * Ctint);
                Spectrum Cspec0 =
                    specular * spec_f0 * (1 - metallic) * spec_color + metallic * base_color;
                Real spec_weight = (1 - specular_transmission * (1 - metallic));

                Spectrum F = schlick_fresnel(Cspec0, h_dot_out);
                Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
                Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
                Real G_out = smith_masking_gtr2(to_local(frame, dir_out), alpha_x, alpha_y);
                Real G = G_in * G_out;
                metallic_contrib = spec_weight * F * D * G / (4 * n_dot_in);
            }

            // Clearcoat component
            if (n_dot_in > 0 && n_dot_h > 0) {
                Real Fc = schlick_fresnel(Real(0.04), h_dot_out);
                // Generalized Trowbridge-Reitz distribution
                Real Dc = GTR1(n_dot_h, alpha_c);
                // SmithG with fixed alpha
                Real Gc_in = smith_masking_gtr2(to_local(frame, dir_in), Real(0.25), Real(0.25));
                Real Gc_out = smith_masking_gtr2(to_local(frame, dir_out), Real(0.25), Real(0.25));
                Real Gc = Gc_in * Gc_out;

                clearcoat_contrib = clearcoat * Fc * Dc * Gc / (4 * n_dot_in);
            }
        }

        // For glass, lights bounce at both sides of the surface.

        // Glass component
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Real Fg = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
        Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
        Real G_out = smith_masking_gtr2(to_local(frame, dir_out), alpha_x, alpha_y);
        Real G = G_in * G_out;
        Spectrum glass_contrib =
            base_color * ((1 - metallic) * specular_transmission * (Fg * D * G) / (4 * fabs(n_dot_in)));

        return diffuse_contrib + sheen_contrib + metallic_contrib + glass_contrib + clearcoat_contrib;
    } else {
        // Only the glass component for refraction
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Vector3 half_vector = normalize(dir_in + dir_out * eta);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;

        Real Fg = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
        Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
        Real G_out = smith_masking_gtr2(to_local(frame, dir_out), alpha_x, alpha_y);
        Real G = G_in * G_out;

        // Burley propose to take the square root of the base color to preserve albedo
        return sqrt(base_color) * ((1 - metallic) * specular_transmission *
            (eta_factor * (1 - Fg) * D * G * eta * eta * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!

    // Fetch the texture values for later use
    Real specular_transmission =
        eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = Real(0.25) * eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);
    Real alpha_c = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Real diffuse_weight = (1 - metallic) * (1 - specular_transmission);
    Real metallic_weight = (1 - specular_transmission * (1 - metallic));
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real clearcoat_weight = clearcoat;
    Real total_weight = diffuse_weight + metallic_weight + glass_weight + clearcoat_weight;
    Real diffuse_prob = diffuse_weight / total_weight;
    Real metallic_prob = metallic_weight / total_weight;
    Real glass_prob = glass_weight / total_weight;
    Real clearcoat_prob = clearcoat_weight / total_weight;

    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // Our incoming ray is coming from inside,
        // so the probability of sampling the glass lobe is 1 if glass_prob is not 0.
        diffuse_prob = 0;
        metallic_prob = 0;
        clearcoat_prob = 0;
        if (glass_prob > 0) {
            glass_prob = 1;
        }
    }

    if (reflect) {
        // For metallic: visible normal sampling -> D * G_in
        Vector3 half_vector = normalize(dir_in + dir_out);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real n_dot_in = dot(frame.n, dir_in);
        Real n_dot_h = dot(half_vector, vertex.shading_frame.n);
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);

        // For diffuse, metallic, and clearcoat, the light bounces
        // only at the top of the surface.
        if (dot(vertex.geometric_normal, dir_in) >= 0 &&
                dot(vertex.geometric_normal, dir_out) >= 0) {
            diffuse_prob *= fmax(dot(frame.n, dir_out), Real(0)) / c_PI;

            if (n_dot_in > 0) {
                Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
                Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
                metallic_prob *= (D * G_in / (4 * n_dot_in));
            } else {
                metallic_prob = 0;
            }

            // For clearcoat: D importance sampling
            if (n_dot_h > 0 && h_dot_out > 0) {
                Real Dc = GTR1(n_dot_h, alpha_c);
                clearcoat_prob *= (Dc * n_dot_h / (4 * h_dot_out));
            } else {
                clearcoat_prob = 0;
            }
        }

        // For glass: F * visible normal
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Real Fg = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
        Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
        glass_prob *= (Fg * D * G_in / (4 * fabs(n_dot_in)));
    } else {
        // Only glass component for refraction
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Vector3 half_vector = normalize(dir_in + dir_out * eta);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);
        Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
        Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
        Real Fg = fresnel_dielectric(h_dot_in, eta);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        glass_prob *= (1 - Fg) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }

    return diffuse_prob + metallic_prob + glass_prob + clearcoat_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!

    // Fetch the texture values for later use
    Real specular_transmission =
        eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = Real(0.25) * eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);
    Real alpha_c = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Real diffuse_weight = (1 - metallic) * (1 - specular_transmission);
    Real metallic_weight = (1 - specular_transmission * (1 - metallic));
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real clearcoat_weight = clearcoat;

    // Two cases: 1) if we are coming from "outside" the surface, 
    // sample all lobes
    if (dot(vertex.geometric_normal, dir_in) >= 0) {
        Real total_weight = diffuse_weight + metallic_weight + glass_weight + clearcoat_weight;
        Real diffuse_prob = diffuse_weight / total_weight;
        Real metallic_prob = metallic_weight / total_weight;
        Real glass_prob = glass_weight / total_weight;
        // Real clearcoat_prob = clearcoat_weight / total_weight;
        if (rnd_param_w <= diffuse_prob) {
            return BSDFSampleRecord{
                to_world(vertex.shading_frame, sample_cos_hemisphere(rnd_param_uv)),
                Real(0) /* eta */, Real(1) /* roughness */};
        } else if (rnd_param_w <= (diffuse_prob + metallic_prob)) { // metallic
            // Visible normal sampling

            // Convert the incoming direction to local coordinates
            Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
            Vector3 local_micro_normal =
                sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

            // Transform the micro normal to world space
            Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
            // Reflect over the world space normal
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected, Real(0) /* eta */, roughness /* roughness */};
        } else if (rnd_param_w <= (diffuse_prob + metallic_prob + glass_prob)) { // glass
            if (glass_prob <= 0) {
                // Just to be safe numerically.
                return {};
            }
            // Visible normal sampling

            // Convert the incoming direction to local coordinates
            Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
            Vector3 local_micro_normal =
                sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

            // Transform the micro normal to world space
            Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
            // Flip half-vector if it's below surface
            if (dot(half_vector, frame.n) < 0) {
                half_vector = -half_vector;
            }

            // Now we need to decide whether to reflect or refract.
            // We do this using the Fresnel term.
            Real h_dot_in = dot(half_vector, dir_in);
            Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
            Real F = fresnel_dielectric(h_dot_in, eta);
            // rescale rnd_param_w from
            // (diffuse_prob + metallic_prob, diffuse_prob + metallic_prob + glass_prob]
            // to
            // (0, 1]
            Real u = (rnd_param_w - (diffuse_prob + metallic_prob)) / glass_prob;
            if (u <= F) {
                // Reflect over the world space normal
                Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
                return BSDFSampleRecord{
                    reflected, Real(0) /* eta */, roughness};
            } else {
                // Refraction
                Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
                if (h_dot_out_sq <= 0) {
                    return {};
                }
                // flip half_vector if needed
                if (h_dot_in < 0) {
                    half_vector = -half_vector;
                }
                Real h_dot_out= sqrt(h_dot_out_sq);
                Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
                return BSDFSampleRecord{refracted, eta, roughness};
            }
        } else { // clearcoat
            // Only importance sampling D

            // Appendix B.2 Burley's note
            Real alpha2 = alpha_c * alpha_c;
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
                reflected, Real(0) /* eta */, sqrt(alpha_c) /* roughness */
            };
        }
    } else {
        // 2) otherwise, only consider the glass lobes.

        // Convert the incoming direction to local coordinates
        Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

        // Transform the micro normal to world space
        Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real h_dot_in = dot(half_vector, dir_in);
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real u = rnd_param_w;
        if (u <= F) {
            // Reflect over the world space normal
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected, Real(0) /* eta */, roughness /* roughness */};
        } else {
            // Refraction
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
