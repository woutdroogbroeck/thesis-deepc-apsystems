classdef InsulinPump < handle
    properties
        param
        basal
    end
    
    methods
        function obj = InsulinPump(param, basal)
            obj.param = param;
            obj.basal = basal;
        end
    
        function insulin = discrete_insulin(obj, amount)

            insulin = amount;
            insulin = round(insulin / obj.param.inc_basal) * obj.param.inc_basal;
            insulin = min(insulin, obj.param.max_basal);
            insulin = max(insulin, obj.param.min_basal);

            
        end

    end
end
