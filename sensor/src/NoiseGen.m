classdef NoiseGen < handle
    properties
        param
        seed
    end
    properties
        noise_seq
        noise_init
        PRECOMPUTE = 10 
        MDL_SAMPLE_TIME = 15  
        e
    end
    methods
        function obj = NoiseGen(param, seed)
            obj.param = param;
            obj.seed = seed;
            rng(seed);
            obj.noise_seq = [];
            obj.noise_init = obj.generate_next_number;
            obj.e = randn;
        end

        function [noise, obj] = generate_noise(obj)
            if isempty(obj.noise_seq)
                obj.noise_seq = obj.get_noise_seq;
            end
            noise = obj.noise_seq(1);
            obj.noise_seq(1) = [];
        end

        function [noise_seq, obj] = get_noise_seq(obj)
            noise15 = [obj.noise_init];
            
            for i = 1:obj.PRECOMPUTE
                noise15 = [noise15, obj.generate_next_number];
            end
            
            obj.noise_init = noise15(end);
            
            t15 = (0:length(noise15)-1) * obj.MDL_SAMPLE_TIME;
            
            nsample = floor(obj.PRECOMPUTE * obj.MDL_SAMPLE_TIME / obj.param.sample_time)+1;
            t = (0:nsample-1) * obj.param.sample_time;
            t = t(t <= max(t15));
            noise = interp1(t15, noise15, t, 'cubic');
            
            noise_seq = noise(2:end);
        end

        function eps = generate_next_number(obj)
            obj.e = obj.param.PACF * (obj.e + randn);
            eps = johnson_transform(obj.param.xi, obj.param.lambda, ...
                obj.param.gamma, obj.param.delta, obj.e);
        end
    end

end