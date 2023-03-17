Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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

    // Disney BSDF's diffuse component contains two parts:
    // 1) a "base" diffuse model based on a modified Fresnel response
    // 2) a "subsurface" diffuse model based on the Lommel-Seeliger law (modified to make things brighter) 
    //    (brought to graphics by Hanrahan & Kruger)

    // Fetch the textures value for later use
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Common variables
    Real n_dot_out = dot(frame.n, dir_out);
    if (n_dot_out <= 0) {
        return make_zero_spectrum();
    }
    Real n_dot_in = dot(frame.n, dir_in);
    Real schlick_n_dot_out = pow(1 - n_dot_out, Real(5));
    Real schlick_n_dot_in  = pow(1 - n_dot_in, Real(5));
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_out = dot(half_vector, dir_out);

    // The base diffuse model
    Real Fd90 = Real(0.5) + 2 * roughness * h_dot_out * h_dot_out;
    Real base_diffuse = (1 + (Fd90 - 1) * schlick_n_dot_out) *
                        (1 + (Fd90 - 1) * schlick_n_dot_in);

    // The subsurface model
    // Disney's hack to increase the response at grazing angle
    Real Fss90 = h_dot_out * h_dot_out * roughness;
    Real Fss = (1 + (Fss90 - 1) * schlick_n_dot_out) *
               (1 + (Fss90 - 1) * schlick_n_dot_in);
    // Lommel-Seeliger law (modified/rescaled)
    Real ss = Real(1.25) * (Fss * (1 / (n_dot_out + n_dot_in) - Real(0.5)) + Real(0.5));

    return base_color * ((base_diffuse * (1 - subsurface) + ss * subsurface) / c_PI) * n_dot_out;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(vertex.shading_frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
