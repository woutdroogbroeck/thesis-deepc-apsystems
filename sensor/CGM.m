classdef CGM < handle
    properties
        param
        Ts
        noise_power
        seed     
        last_Gm
        N
    end
    properties
        noise
        lastN_measurements = []
        e
    end
    methods
        function obj = CGM(param, noise_power, seed, N)
            obj.param = param;
            obj.Ts = param.sample_time;
            obj.seed = seed;
            obj.noise = NoiseGen(param, seed);
            obj.noise_power = noise_power;
            obj.N = N;
            rng(seed);
            obj.e = rand;
        end

        function [Gm, obj] = measure(obj, patient)
            if (mod(patient.t, obj.Ts) == 0)
                G = patient.G;
                Gm = G + obj.noise_power*obj.noise.generate_noise;
                Gm = max(Gm, obj.param.min);
                Gm = min(Gm, obj.param.max);
                obj.last_Gm = Gm;
            else
                Gm = obj.last_Gm;
            end  
            obj.lastN_measurements = [obj.lastN_measurements, Gm];
            if length(obj.lastN_measurements) > obj.N
                obj.lastN_measurements(1) = [];
            end
        end
        
        function obj = set_seed(obj, seed)
            obj.seed = seed;
        end

        
        function obj = reset(obj)
            obj.last_Gm = 0;
            obj.noise = NoiseGen(obj.param, obj.seed);
            obj.lastN_measurements = [];
        end
    end
end