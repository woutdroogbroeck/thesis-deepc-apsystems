function patient = load_patient(model, Ts, patient_id, linear)
%load_patient load a patient instance into the workspace
%   model : "minimal" or "uvapadova"
%   Ts    : sample time / discretization step
%   id    : patient id -> minimal: 1-3 | uvapadova: 1-30
%   returns a patient object

    if model == "minimal"
        all_param = load('patient/parameters/minimal_patients.mat').minimal_patients;
        param = all_param(patient_id,:);

    elseif model == "uvapadova"
      
        all_param = load('patient/parameters/maximal_patients.mat').maximal_patients;
        param = all_param(patient_id,:);
    end
    
    patient = Patient(param,Ts,model,linear); 
end