clear
script2protect = {"closure_depth6_prior_model1227.m";
    "closure_temps.m";
    "closure_temps_linear.m";
    "inversion_clc2.m";
    "isotherms_clc0909.m";
    "kvalue.m";
    "linexline.m";
    "pred_staris.m";
    "temp_mean.m";
    "tridag.m"};
for i=1:size(script2protect,1)
    pcode(script2protect{i,1});
    delete(script2protect{i,1})
end
