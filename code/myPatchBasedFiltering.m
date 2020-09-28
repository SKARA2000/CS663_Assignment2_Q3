%% Optimizing function
function [OptimalRMSD] = myPatchBasedFiltering(path, format, name, tune)
    if(format=="mat")
        load(path, "image")
        image = im2double(image);
    elseif(format=="png")
        image = im2double(imread(path,'png'));
    elseif(format=="im")
        image = im2double(path);
    end
    [m3,n3] = size(image);
    rng('default');
    Noisyimage = image + (randn(m3,n3)*0.05);
    rmse=10;
    optimal_sd1=1.3;
    if(name=="barbara")
        optimal_sd2=0.044;
        savepath="../images/barbara";
    elseif(name=="grass")
        optimal_sd2=0.04;
        savepath="../images/grass";
    elseif(name=="honeycomb")
        optimal_sd2=0.044;
        savepath="../images/honeycomb";
    end
    if(tune=='y')
        if(name=="barbara")
            for sd2=0.04:0.001:0.05
                   sample_rmsd = myPatchBasedFilteringImg(Noisyimage, image, image,optimal_sd1, sd2, 'n', savepath);
                   if(sample_rmsd<rmse)
                       %%optimal_sd1 = sd1;
                       optimal_sd2 = sd2;
                       rmse = sample_rmsd;
                   end
            end
        elseif(name=="grass")
            for sd2=0.03:0.005:0.06
                   sample_rmsd = myPatchBasedFilteringImg(Noisyimage, image, optimal_sd1, sd2, 'n',savepath);
                   if(sample_rmsd<rmse)
                       %%optimal_sd1 = sd1;
                       optimal_sd2 = sd2;
                       rmse = sample_rmsd;
                   end
            end
        elseif(name=="honeycomb")
            for sd2=0.03:0.002:0.06
                   sample_rmsd = myPatchBasedFilteringImg(Noisyimage, image, optimal_sd1, sd2, 'n',savepath);
                   if(sample_rmsd<rmse)
                       %%optimal_sd1 = sd1;
                       optimal_sd2 = sd2;
                       rmse = sample_rmsd;
                   end
            end
        end
    end
    sample_patch_kernel = constructGaussian(optimal_sd1, 9);
    figure
    imagesc(sample_patch_kernel);
    colormap gray;
    title(["Patch Gaussian Kernel for ",name]);
    colorbar;
    savefig(savepath+"_patch_kernel.fig");
    image1=myPatchBasedFilteringImg(Noisyimage, image, optimal_sd1, optimal_sd2, 'y',savepath+1+".fig");
    if(tune=='n')
        image2=myPatchBasedFilteringImg(Noisyimage, image, 0.9*optimal_sd1, optimal_sd2, 'y',savepath+2+".fig");
        image3=myPatchBasedFilteringImg(Noisyimage, image, 1.1*optimal_sd1, optimal_sd2, 'y',savepath+3+".fig");
        image4=myPatchBasedFilteringImg(Noisyimage, image, optimal_sd1, 0.9*optimal_sd2, 'y',savepath+4+".fig");
        image5=myPatchBasedFilteringImg(Noisyimage, image, optimal_sd1, 1.1*optimal_sd2, 'y',savepath+5+".fig");
    end
    OptimalRMSD = image1;
    fprintf(name);
    fprintf(":\n")
    fprintf("Optimal RMSE=%d \n",image1);
    if(tune=='n')
        fprintf("RMSE for 0.9*sd_space=%d \n",image2);
        fprintf("RMSE for 1.1*sd_space=%d \n",image3);
        fprintf("RMSE for 0.9*sd_intensity=%d \n",image4);
        fprintf("RMSE for 1.1*sd_intensity=%d \n",image5);
    end
end


%% Function to implement Patch Based Filtering
function [RMSD] = myPatchBasedFilteringImg(imageOrig, imageProper, sd_patch, sd_kernel, show, savepath)
    %{
    sd_patch = 1.3;
    sd_kernel = 0.01;
    %}
    %%till now optimal 0.04
    window_size = 25;
    patch_kernel_size = 9;
    [m,n] = size(imageOrig);
    rng('default');
    ImageBlur = imageOrig;
    patch_kernel = constructGaussian(sd_patch,patch_kernel_size);
    s_half = int16((window_size+patch_kernel_size)/2);
    p_half = int16((patch_kernel_size-1)/2);
    k_half = int16((window_size-1)/2);
    Numerator_kernel=0;
    Denominator_kernel=0;
    for i=1:m
        for j=1:n
            if((i>=s_half) && (i<=(m-s_half+1)) && (j>=s_half) && (j<=(n-s_half+1)))
                window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                for a=(i-k_half):(i+k_half)
                    for b=(j-k_half):(j+k_half)
                        patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                        term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                        sum_term = sum(term_patch,'all'); 
                        Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                        Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                    end
                end
            elseif((i<s_half) && (i>=(p_half+1)))
                if((j<s_half) && (j>=(p_half+1)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(p_half+1):(i+k_half)
                        for b=(p_half+1):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>(n-s_half+1)) && (j<=(n-p_half)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(p_half+1):(i+k_half)
                        for b=(j-k_half):(n-p_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>=s_half) && (j<=(n-s_half+1)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(p_half+1):(i+k_half)
                        for b=(j-k_half):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>(n-p_half)) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                 elseif((j<=p_half) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                end
            elseif((i>(m-s_half+1)) && (i<=(m-p_half)))
                if((j<s_half) && (j>=(p_half+1)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(i-k_half):(m-p_half)
                        for b=(p_half+1):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>(n-s_half+1)) && (j<=(n-p_half)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(i-k_half):(m-p_half)
                        for b=(j-k_half):(n-p_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>=s_half) && (j<=(n-s_half+1)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(i-k_half):(m-p_half)
                        for b=(j-k_half):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>(n-p_half)) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                 elseif((j<=p_half) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                end
            elseif((i>=s_half) && (i<=(m-s_half+1)))
                if((j<s_half) && (j>=(p_half+1)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(i-k_half):(i+k_half)
                        for b=(p_half+1):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>(n-s_half+1)) && (j<=(n-p_half)))
                    window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):(j+p_half));
                    for a=(i-k_half):(i+k_half)
                        for b=(j-k_half):(n-p_half)
                            patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*patch_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end 
                elseif((j>(n-p_half)) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),(j-p_half):end);
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,1:n2);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=(j-k_half):j
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b-p_half):(n2+b-p_half-1));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                 elseif((j<=p_half) && (i>p_half) && (i<=(m-p_half)))
                    if((i>p_half) && (i<=(k_half+p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(p_half+1):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    elseif((i<=(m-p_half)) && (i>(m-k_half-p_half)))
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(m-p_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    else
                        window_orig = imageOrig((i-p_half):(i+p_half),1:(j+p_half));
                        [m2,n2] = size(window_orig);
                        work_kernel = patch_kernel(1:end,(patch_kernel_size-n2+1):end);
                        work_kernel = work_kernel./(sum(work_kernel,'all'));
                        for a=(i-k_half):(i+k_half)
                            for b=j:(j+k_half)
                                patch_sample = imageOrig((a-p_half):(a+p_half),(b+p_half-n2+1):(b+p_half));
                                term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                                sum_term = sum(term_patch,'all'); 
                                Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                                Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                            end
                        end
                    end
                end
            elseif((i<=p_half))
                if(j<=p_half)
                    window_orig = imageOrig(1:(i+p_half),1:(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel((patch_kernel_size-m2+1):end,(patch_kernel_size-n2+1):end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=i:(i+k_half)
                        for b=j:(j+k_half)
                            patch_sample = imageOrig((a+p_half-m2+1):(a+p_half),(b+p_half-n2+1):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif(j>(n-p_half))
                    window_orig = imageOrig(1:(i+p_half),(j-p_half):end);
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel((patch_kernel_size-m2+1):end,1:n2);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=i:(i+k_half)
                        for b=(j-k_half):j
                            patch_sample = imageOrig((a+p_half-m2+1):(a+p_half),(b-p_half):(n2+b-p_half-1));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>p_half) && (j<=(k_half+p_half)))
                    window_orig = imageOrig(1:(i+p_half),(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel((patch_kernel_size-m2+1):end,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=i:(i+k_half)
                        for b=(p_half+1):(j+k_half)
                            patch_sample = imageOrig((a+p_half-m2+1):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j<=(n-p_half)) && (j>(n-k_half-p_half)))
                    window_orig = imageOrig(1:(i+p_half),(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel((patch_kernel_size-m2+1):end,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=i:(i+k_half)
                        for b=(j-k_half):(n-p_half)
                            patch_sample = imageOrig((a+p_half-m2+1):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                else
                    window_orig = imageOrig(1:(i+p_half),(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel((patch_kernel_size-m2+1):end,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=i:(i+k_half)
                        for b=(j-k_half):(j+k_half)
                            patch_sample = imageOrig((a+p_half-m2+1):(a+p_half),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                end
            elseif(i>(m-p_half))
                if(j<=p_half)
                    window_orig = imageOrig((i-p_half):end,1:(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel(1:m2,(patch_kernel_size-n2+1):end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=(i-k_half):i
                        for b=j:(j+k_half)
                            patch_sample = imageOrig((a-p_half):(m2+a-p_half-1),(b+p_half-n2+1):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif(j>(n-p_half))
                    window_orig = imageOrig((i-p_half):end,(j-p_half):end);
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel(1:m2,1:n2);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=(i-k_half):i
                        for b=(j-k_half):j
                            patch_sample = imageOrig((a-p_half):(m2+a-p_half-1),(b-p_half):(n2+b-p_half-1));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j>p_half) && (j<=(k_half+p_half)))
                    window_orig = imageOrig((i-p_half):end,(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel(1:m2,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=(i-k_half):i
                        for b=(p_half+1):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(m2+a-p_half-1),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                elseif((j<=(n-p_half)) && (j>(n-k_half-p_half)))
                    window_orig = imageOrig((i-p_half):end,(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel(1:m2,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=(i-k_half):i
                        for b=(j-k_half):(n-p_half)
                            patch_sample = imageOrig((a-p_half):(m2+a-p_half-1),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                else
                    window_orig = imageOrig((i-p_half):end,(j-p_half):(j+p_half));
                    [m2,n2] = size(window_orig);
                    work_kernel = patch_kernel(1:m2,1:end);
                    work_kernel = work_kernel./(sum(work_kernel,'all'));
                    for a=(i-k_half):i
                        for b=(j-k_half):(j+k_half)
                            patch_sample = imageOrig((a-p_half):(m2+a-p_half-1),(b-p_half):(b+p_half));
                            term_patch = (window_orig - patch_sample).*(window_orig - patch_sample).*work_kernel;
                            sum_term = sum(term_patch,'all'); 
                            Numerator_kernel = Numerator_kernel + (gaussian_distance(sum_term,sd_kernel)*imageOrig(a,b));
                            Denominator_kernel = Denominator_kernel + gaussian_distance(sum_term,sd_kernel);
                        end
                    end
                end
            end
            ImageBlur(i,j) = Numerator_kernel/Denominator_kernel;
            Numerator_kernel=0;
            Denominator_kernel=0;
        end
    end
    RMSD = sqrt(sum((imageProper-ImageBlur).*(imageProper-ImageBlur),'all')/(m*n));
    if(show=='y')
        figure
        subplot(1,3,1), imshow(imageProper);
        colorbar
        axis on
        title("Original Image")
        subplot(1,3,2), imshow(imageOrig);
        colorbar
        axis on
        title("Corrupted Image")
        subplot(1,3,3), imshow(ImageBlur);
        colorbar
        axis on
        title({'Filtered Image'; ['\sigma_{patch}=',num2str(sd_patch),', \sigma_{kernel}=',num2str(sd_kernel)]});
        savefig(savepath);
    end
end
%%
function [Gaussian] = constructGaussian(sd,size)
    ctr = (size+1)/2;
    Gaussian = zeros(size, size);
    for i = 1:size
        for j = 1:size
            Gaussian(i,j) = gaussian((i-ctr),(j-ctr));
        end
    end
    Gaussian = Gaussian.*sd;
    Gaussian = Gaussian./sum(Gaussian,'all');
end

function val = gaussian(x,y)
    val = (1/sqrt(2*pi))*exp(-((x^2)+(y^2))/(2));
end

function gauss_Value = gaussian_distance(x,sd_gauss_general)
    gauss_Value = (1/sqrt(2*pi)*sd_gauss_general)*exp(-(x)/(2*(sd_gauss_general^2)));
end

        