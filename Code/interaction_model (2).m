%%Cargar modelo de astrocito
astrocito = readCbModel('Astrocyte_Mendoza2022.xml');
%%%Cargar modelo de neurona
loaded_model = load('modelo_neurona_control_estequiometric_biomass.mat');
neurona = loaded_model.modelo_control;
%%% normalizar los nomres de metabolitos y reacciones

%% 1. Prefijar modelos
astrocyteModel = prefixModelFields(modelo_control_astrocito, 'Astro_');
neuronModel    = prefixModelFields(neurona, 'Neuron_');

%% 2. Concatenar componentes del modelo
mergedModel.rxns  = [astrocyteModel.rxns; neuronModel.rxns];
mergedModel.mets  = [astrocyteModel.mets; neuronModel.mets];
mergedModel.genes = [astrocyteModel.genes; neuronModel.genes];

% Concatenar matriz estequiométrica con ceros cruzados
mergedModel.S = [astrocyteModel.S, sparse(size(astrocyteModel.S,1), size(neuronModel.S,2));
                 sparse(size(neuronModel.S,1), size(astrocyteModel.S,2)), neuronModel.S];

% Concatenar vectores de límites y función objetivo
mergedModel.lb = [astrocyteModel.lb; neuronModel.lb];
mergedModel.ub = [astrocyteModel.ub; neuronModel.ub];
mergedModel.c  = [astrocyteModel.c;  neuronModel.c];

% Concatenar fórmulas químicas si están disponibles
if isfield(astrocyteModel, 'metFormulas') && isfield(neuronModel, 'metFormulas')
    mergedModel.metFormulas = [astrocyteModel.metFormulas; neuronModel.metFormulas];
end

%% 3. Eliminar duplicados
[~, uniqueMetIdx] = unique(mergedModel.mets, 'stable');
mergedModel.mets = mergedModel.mets(uniqueMetIdx);

[~, uniqueRxnIdx] = unique(mergedModel.rxns, 'stable');
mergedModel.rxns = mergedModel.rxns(uniqueRxnIdx);

% Filtrar matriz S según nuevos índices
mergedModel.S = mergedModel.S(uniqueMetIdx, uniqueRxnIdx);

% Ajustar lb, ub, c a las nuevas reacciones
mergedModel.lb = mergedModel.lb(uniqueRxnIdx);
mergedModel.ub = mergedModel.ub(uniqueRxnIdx);
mergedModel.c  = mergedModel.c(uniqueRxnIdx);

%% 4. Asegurar dimensiones consistentes
if length(mergedModel.lb) ~= length(mergedModel.rxns)
    mergedModel.lb = mergedModel.lb(1:length(mergedModel.rxns));
end
if length(mergedModel.ub) ~= length(mergedModel.rxns)
    mergedModel.ub = mergedModel.ub(1:length(mergedModel.rxns));
end
if length(mergedModel.c) ~= length(mergedModel.rxns)
    mergedModel.c = zeros(length(mergedModel.rxns), 1);
end
if size(mergedModel.S, 1) ~= length(mergedModel.mets) || ...
   size(mergedModel.S, 2) ~= length(mergedModel.rxns)
    error('La matriz S no coincide con el número de metabolitos o reacciones.');
end

%% 5. Optimizar modelo
solution = optimizeCbModel(mergedModel);

if solution.stat == 1
    disp('Optimización exitosa.');
    disp(['Valor de la función objetivo: ', num2str(solution.f)]);
else
    disp('El modelo no pudo optimizarse.');
end

%% 6. Verificar balance de masa y carga (opcional)
if isfield(mergedModel, 'metFormulas') && ~isempty(mergedModel.metFormulas)
    try
        mergedModel = computeMetFormulaCharge(mergedModel);
        [balancedRxns, unbalancedRxns] = checkMassChargeBalance(mergedModel);

        if ~isempty(unbalancedRxns)
            disp('Reacciones no balanceadas:');
            disp(unbalancedRxns);
        else
            disp('Todas las reacciones están balanceadas.');
        end
    catch ME
        warning('No se pudo verificar el balance de masa y carga: %s', ME.message);
    end
else
    disp('Advertencia: No hay fórmulas químicas disponibles para verificar el balance.');
end

%% 7. Guardar modelo
%save('mergedModel_control.mat', 'mergedModel');

% Mostrar resumen
fprintf('\nResumen del modelo fusionado:\n');
fprintf('Reacciones: %d\n', length(mergedModel.rxns));
fprintf('Metabolitos: %d\n', length(mergedModel.mets));
fprintf('Genes: %d\n', length(mergedModel.genes));
fprintf('Tamaño de S: %d x %d\n', size(mergedModel.S,1), size(mergedModel.S,2));

% Obtener índices de las funciones objetivo individuales
idxAstroBiomass  = find(strcmp(mergedModel.rxns, 'Astro_biomass_maintenance'));
idxNeuronBiomass = find(strcmp(mergedModel.rxns, 'Neuron_HMR_0308'));

% Reiniciar el vector de la función objetivo
mergedModel.c = zeros(length(mergedModel.rxns), 1);

% Definir la función objetivo comunitaria (suma de ambas)
mergedModel.c(idxAstroBiomass)  = 1;
mergedModel.c(idxNeuronBiomass) = 1;

% Optimizar el modelo con función objetivo comunitaria
solution = optimizeCbModel(mergedModel);

% Mostrar los resultados
disp(['Valor de la biomasa comunitaria: ', num2str(solution.f)]);
disp(['Flujo Astro_biomass_maintenance: ', num2str(solution.x(idxAstroBiomass))]);
disp(['Flujo Neuron_HMR_0308: ', num2str(solution.x(idxNeuronBiomass))]);


mergeModel= load('mergedModel_control.mat')

%%%%%%
%El modelo que tienes en este momento es una fusión de los modelos de astrocito y neurona, 
% pero para que sea considerado un modelo de interacción entre astrocito y neurona, debes asegurarte 
% de que los dos tipos celulares estén intercambiando metabolitos clave de manera explícita, lo que 
% representa una interacción metabólica entre ambas células.

% Obtener los nombres de los metabolitos
metaboliteNames = mergedModel.mets;

% Extraer los compartimentos de los nombres de los metabolitos
compartments = regexp(metaboliteNames, '\[.*?\]', 'match');
compartments = unique([compartments{:}]);

% Mostrar los compartimentos únicos en el modelo
disp('Compartimentos en el modelo fusionado:');
disp(compartments);


% Verificar si el compartimento extracelular está presente
if any(strcmp(compartments, '[e]'))
    disp('El compartimento extracelular está presente.');
else
    disp('El compartimento extracelular no está presente. Debe añadirse.');
end

% Buscar metabolitos en el compartimento extracelular '[e]'
extracellularMets = mergedModel.mets(contains(mergedModel.mets, '[e]'));

% Mostrar los metabolitos en el compartimento extracelular
disp('Metabolitos en el compartimento extracelular:');
disp(extracellularMets);

% Verificar reacciones de intercambio con el compartimento extracelular
exchangeReactions = mergedModel.rxns(contains(mergedModel.rxns, 'EX_'));

% Mostrar las reacciones de intercambio
disp('Reacciones de intercambio en el compartimento extracelular:');
disp(exchangeReactions);

%el siguiente paso es asegurarte de que estas reacciones están bien configuradas 
% para permitir el flujo de metabolitos entre el compartimento extracelular y los compartimentos de los astrocitos y las neuronas.

% Verificar los límites inferiores y superiores de las reacciones de intercambio
for i = 1:length(exchangeReactions)
    rxnID = find(strcmp(mergedModel.rxns, exchangeReactions{i}));
    disp(['Reacción: ', exchangeReactions{i}, ' lb: ', num2str(mergedModel.lb(rxnID)), ', ub: ', num2str(mergedModel.ub(rxnID))]);
end

%%%%%___________________________________________
% Prefijo para los metabolitos de astrocito en varios compartimentos
astroCompartments = {'[c]', '[m]', '[g]', '[r]', '[e]'};  % Especificar aquí los compartimentos relevantes para el astrocito
for i = 1:length(astroCompartments)
    astrocyteMetIdx = find(contains(mergedModel.mets, astroCompartments{i}));
    mergedModel.mets(astrocyteMetIdx) = strcat('Astro_', mergedModel.mets(astrocyteMetIdx));
end


% Prefijo para los metabolitos de neurona en varios compartimentos
neuronCompartments = {'[c]', '[m]', '[g]', '[r]'};  % Compartimentos relevantes para la neurona
for i = 1:length(neuronCompartments)
    neuronMetIdx = find(contains(mergedModel.mets, neuronCompartments{i}) & ~contains(mergedModel.mets, 'Astro_'));  % Evitar doble prefijo
    mergedModel.mets(neuronMetIdx) = strcat('Neuron_', mergedModel.mets(neuronMetIdx));
end

% Mostrar los metabolitos que pertenecen al astrocito
disp('Metabolitos del astrocito:');
disp(mergedModel.mets(contains(mergedModel.mets, 'Astro_')));

% Mostrar los metabolitos que pertenecen a la neurona
disp('Metabolitos de la neurona:');
disp(mergedModel.mets(contains(mergedModel.mets, 'Neuron_')));


% Verificar las reacciones de intercambio en el compartimento extracelular
exchangeReactions = mergedModel.rxns(contains(mergedModel.rxns, 'EX_') & contains(mergedModel.rxns, '[e]'));
disp('Reacciones de intercambio en el compartimento extracelular:');
disp(exchangeReactions);

% Reacciones de interacción astrocito-neurona
reacciones_interaccion = {
    'Astro_EX_lac_L[e]', 'Neuron_EX_lac_L[e]'
    'Astro_EX_glu_L[e]', 'Neuron_EX_glu_L[e]'
    'Astro_EX_ala_D[e]', 'EX_ala_D[e]'
    'Astro_EX_gg4abut[e]', 'EX_gg4abut[e]'
    'Astro_EX_asp_D[e]', 'EX_asp_D[e]'
};

disp('--- Flujos de metabolitos clave ---');
for i = 1:size(reacciones_interaccion,1)
    for j = 1:2
        rxn = reacciones_interaccion{i,j};
        if any(strcmp(mergedModel.rxns, rxn))
            flujo = solution.x(strcmp(mergedModel.rxns, rxn));
            disp([rxn, ': ', num2str(flujo)]);
        else
            disp([rxn, ' no está presente en el modelo.']);
        end
    end
end

% Lista de nuevas reacciones de importación astrocítica
importList = {
    'Astro_lac_L_in',     {'lac_L[e]', 'lac_L[c]'},       [-1, 1];
    'Astro_glu_L_in',     {'glu_L[e]', 'glu_L[c]'},       [-1, 1];
    'Astro_ala_D_in',     {'ala_D[e]', 'ala_D[c]'},       [-1, 1];
    'Astro_gg4abut_in',   {'gg4abut[e]', 'gg4abut[c]'},   [-1, 1];
    'Astro_asp_D_in',     {'asp_D[e]', 'asp_D[c]'},       [-1, 1];
};

% Crear las reacciones que no estén en el modelo
for i = 1:size(importList, 1)
    rxnID = importList{i,1};
    mets = importList{i,2};
    stoich = importList{i,3};

    if any(strcmp(mergedModel.rxns, rxnID))
        fprintf('✔ Ya existe: %s\n', rxnID);
        continue;
    end

    % Crear columna nueva en S
    mergedModel.S(:,end+1) = sparse(length(mergedModel.mets), 1);
    mergedModel.rxns{end+1} = rxnID;

    for j = 1:length(mets)
        metIndex = find(strcmp(mergedModel.mets, mets{j}));
        if isempty(metIndex)
            warning('⚠️ Metabolito no encontrado: %s\n', mets{j});
        else
            mergedModel.S(metIndex, end) = stoich(j);
        end
    end

    % Limites y función objetivo
    mergedModel.lb(end+1) = 0;
    mergedModel.ub(end+1) = 1000;
    mergedModel.c(end+1)  = 0;

    % Vacío para reglas
    mergedModel.grRules{end+1} = '';
    mergedModel.rules{end+1}   = '';
end

% Optimización del modelo
solution = optimizeCbModel(mergedModel);

% Mostrar resultado general
disp(['Valor de la función objetivo: ', num2str(solution.f)]);

% Mostrar flujos de importación astrocítica
rxnsToCheck = {
    'Astro_lac_L_in', 
    'Astro_glu_L_in', 
    'Astro_ala_D_in', 
    'Astro_gg4abut_in',
    'Astro_asp_D_in'
};

for i = 1:length(rxnsToCheck)
    idx = find(strcmp(mergedModel.rxns, rxnsToCheck{i}));
    if ~isempty(idx)
        fprintf('Flujo %s: %.4f\n', rxnsToCheck{i}, solution.x(idx));
    else
        fprintf('⚠️ Reacción %s no encontrada\n', rxnsToCheck{i});
    end
end

% --- Exportación de Alanina D ---
met_ala_c = 'ala_D[c]';
met_ala_e = 'ala_D[e]';
if all(ismember({met_ala_c, met_ala_e}, mergedModel.mets))
    idx_c = find(strcmp(mergedModel.mets, met_ala_c));
    idx_e = find(strcmp(mergedModel.mets, met_ala_e));
    
    mergedModel.S(:, end+1) = sparse(length(mergedModel.mets), 1);
    mergedModel.S(idx_c, end) = -1;
    mergedModel.S(idx_e, end) = 1;
    
    mergedModel.rxns{end+1} = 'Neuron_EX_ala_D[e]';
    mergedModel.lb(end+1) = 0;
    mergedModel.ub(end+1) = 1000;
    mergedModel.c(end+1) = 0;
    disp('✅ Reacción añadida: Neuron_EX_ala_D[e]');
else
    warning('⚠️ Metabolitos de alanina D no encontrados.');
end

% --- Exportación de Aspartato D ---
met_asp_c = 'asp_D[c]';
met_asp_e = 'asp_D[e]';
if all(ismember({met_asp_c, met_asp_e}, mergedModel.mets))
    idx_c = find(strcmp(mergedModel.mets, met_asp_c));
    idx_e = find(strcmp(mergedModel.mets, met_asp_e));
    
    mergedModel.S(:, end+1) = sparse(length(mergedModel.mets), 1);
    mergedModel.S(idx_c, end) = -1;
    mergedModel.S(idx_e, end) = 1;
    
    mergedModel.rxns{end+1} = 'Neuron_EX_asp_D[e]';
    mergedModel.lb(end+1) = 0;
    mergedModel.ub(end+1) = 1000;
    mergedModel.c(end+1) = 0;
    disp('✅ Reacción añadida: Neuron_EX_asp_D[e]');
else
    warning('⚠️ Metabolitos de aspartato D no encontrados.');
end

% Reoptimizar
solution = optimizeCbModel(mergedModel);
disp(['\n🔁 Valor de la función objetivo: ', num2str(solution.f)]);


% --- Verificación y adición de metabolitos extracelulares e intracelulares clave ---
requiredMets = {
    'lac_L[e]',  'C3H5O3',     -1
    'lac_L[c]',  'C3H5O3',     -1
    'glu_L[e]',  'C5H9NO4',     0
    'glu_L[c]',  'C5H9NO4',     0
    'ala_D[e]',  'C3H7NO2',     0
    'ala_D[c]',  'C3H7NO2',     0
    'gg4abut[e]','C4H9NO2',     0
    'gg4abut[c]','C4H9NO2',     0
    'asp_D[e]',  'C4H7NO4',    -1
    'asp_D[c]',  'C4H7NO4',    -1
};

for i = 1:size(requiredMets,1)
    metName = requiredMets{i,1};
    metFormula = requiredMets{i,2};
    metCharge  = requiredMets{i,3};

    if ~any(strcmp(mergedModel.mets, metName))
        mergedModel.mets{end+1} = metName;

        % Añadir a metFormulas y metCharges si existen los campos
        if isfield(mergedModel, 'metFormulas')
            mergedModel.metFormulas{end+1} = metFormula;
        end
        if isfield(mergedModel, 'metCharges')
            mergedModel.metCharges(end+1) = metCharge;
        end

        % Añadir nueva fila vacía en S
        mergedModel.S(end+1, :) = sparse(1, length(mergedModel.rxns));
        fprintf('✅ Metabolito añadido: %s\n', metName);
    else
        fprintf('✔ Metabolito ya existe: %s\n', metName);
    end
end


%%%%%%%%%%%%%%%_________________________________
% Definir los metabolitos involucrados
metabolites = {'lac_L[e]', 'lac_L[c]'};

% Definir los coeficientes estequiométricos correspondientes
stoichCoeffs = [-1, 1];  % Salida de lac_L[e], entrada de lac_L[c]

% Añadir la nueva reacción manualmente
newReactionID = 'Neuron_lac_L_in';  % Nombre de la nueva reacción
mergedModel.rxns{end+1} = newReactionID;  % Añadir la reacción al campo de reacciones
mergedModel.S(:, end+1) = sparse(length(mergedModel.mets), 1);  % Añadir una columna vacía a la matriz estequiométrica

% Añadir las estequiometrías a la matriz S
for i = 1:length(metabolites)
    metIndex = find(strcmp(mergedModel.mets, metabolites{i}));
    mergedModel.S(metIndex, end) = stoichCoeffs(i);
end

% Definir los límites de la reacción
mergedModel.lb(end+1) = 0;  % Límite inferior (importación)
mergedModel.ub(end+1) = 1000;  % Límite superior (importación)

% Añadir una función objetivo nula para la nueva reacción
mergedModel.c(end+1) = 0;

% Definir los metabolitos involucrados para la importación de glutamato
metabolites_glu_import = {'glu_L[e]', 'glu_L[c]'};

% Definir los coeficientes estequiométricos correspondientes
stoichCoeffs_glu_import = [-1, 1];  % Importación de glu_L[e] hacia glu_L[c]

% Añadir la nueva reacción manualmente
newReactionID_glu_import = 'Astro_glu_L_import';  % Nombre de la nueva reacción
mergedModel.rxns{end+1} = newReactionID_glu_import;  % Añadir la reacción al campo de reacciones
mergedModel.S(:, end+1) = sparse(length(mergedModel.mets), 1);  % Añadir una columna vacía a la matriz estequiométrica

% Añadir las estequiometrías a la matriz S
for i = 1:length(metabolites_glu_import)
    metIndex = find(strcmp(mergedModel.mets, metabolites_glu_import{i}));
    mergedModel.S(metIndex, end) = stoichCoeffs_glu_import(i);
end

% Definir los límites de la reacción
mergedModel.lb(end+1) = 0;  % Límite inferior (importación)
mergedModel.ub(end+1) = 1000;  % Límite superior (importación)

% Añadir una función objetivo nula para la nueva reacción
mergedModel.c(end+1) = 0;

% Verificar si las reacciones de glutamato se añadieron correctamente
gluRxns = mergedModel.rxns(contains(mergedModel.rxns, 'glu_L'));
disp('Reacciones relacionadas con glutamato:');
disp(gluRxns);

% Ajustar los límites de flujo para la exportación de glutamato por las neuronas
mergedModel = changeRxnBounds(mergedModel, 'Neuron_glu_L_export', 10, 'u');  % Exportación de glutamato

% Ajustar los límites de flujo para la importación de glutamato en los astrocitos
mergedModel = changeRxnBounds(mergedModel, 'Astro_glu_L_import', 10, 'u');  % Importación de glutamato


% Paso 1: Definir los metabolitos involucrados
metabolites = {'ala_D[e]', 'ala_D[c]'};  % Alanina extracelular e intracelular

% Paso 2: Definir los coeficientes estequiométricos correspondientes
stoichCoeffs = [-1, 1];  % Salida de ala_D[e], entrada de ala_D[c] (captación por astrocitos)

% Paso 3: Añadir la nueva reacción para la captación de alanina D en astrocitos
newReactionID = 'Astro_ala_D_in';  % Nombre de la nueva reacción
mergedModel.rxns{end+1} = newReactionID;  % Añadir la reacción al campo de reacciones
mergedModel.S(:, end+1) = sparse(length(mergedModel.mets), 1);  % Añadir una columna vacía a la matriz estequiométrica

% Paso 4: Añadir las estequiometrías a la matriz S
for i = 1:length(metabolites)
    metIndex = find(strcmp(mergedModel.mets, metabolites{i}));
    mergedModel.S(metIndex, end) = stoichCoeffs(i);
end

% Paso 5: Definir los límites de la reacción
mergedModel.lb(end+1) = 0;  % Límite inferior (importación)
mergedModel.ub(end+1) = 100;  % Límite superior para importación

% Paso 6: Añadir una función objetivo nula para la nueva reacción
mergedModel.c(end+1) = 0;

% Paso 7: Optimizar el modelo para simular la interacción entre astrocitos y neuronas
solution = optimizeCbModel(mergedModel);

% Paso 8: Mostrar el resultado de la optimización
disp(['Valor de la función objetivo: ', num2str(solution.f)]);

% Verificar el flujo de alanina captada por los astrocitos
alanineFluxAstro = solution.x(find(strcmp(mergedModel.rxns, 'Astro_ala_D_in')));
disp(['Flujo de alanina D en la reacción Astro_ala_D_in (astrocitos): ', num2str(alanineFluxAstro)]);

% Verificar el flujo de alanina exportada por las neuronas
alanineFluxNeuron = solution.x(find(strcmp(mergedModel.rxns, 'EX_ala_D[e]')));
disp(['Flujo de alanina D en la reacción EX_ala_D[e] (neuronas): ', num2str(alanineFluxNeuron)]);


mergedModel = changeObjective(mergedModel, 'biomass_reaction');  % Cambia 'BIOMASS_reaction_ID' por la reacción de biomasa
% Optimizar el modelo para simular la interacción entre astrocitos y neuronas
solution = optimizeCbModel(mergedModel);

%Mostrar el resultado de la optimización
dis% Verificar dimensiones de la matriz S
disp('Dimensiones de la matriz S:');
disp(size(mergedModel.S));

% Verificar los metabolitos involucrados
disp('Metabolitos en el modelo:');
disp(mergedModel.mets);

% Verificar límites y coeficientes
disp('Límites inferiores:');
disp(mergedModel.lb);
disp('Límites superiores:');
disp(mergedModel.ub);
disp('Coeficientes de la función objetivo:');
disp(mergedModel.c);

% Intentar optimizar nuevamente sin la nueva reacción
solution = optimizeCbModel(mergedModel);
disp(['Valor de la función objetivo: ', num2str(solution.f)]);
% Verificar el flujo a través de la reacción de importación de lactato en neuronas
lactateFluxNeuron = solution.x(find(strcmp(mergedModel.rxns, 'Neuron_lac_L_in')));
disp(['Flujo de lactato en la reacción Neuron_lac_L_in: ', num2str(lactateFluxNeuron)]);


% Mostrar el resultado de la optimización
disp(['Valor de la función objetivo: ', num2str(solution.f)]);

% Verificar el flujo a través de la reacción de importación de lactato en neuronas
lactateFluxNeuron = solution.x(find(strcmp(mergedModel.rxns, 'Neuron_lac_L_in')));
disp(['Flujo de lactato en la reacción Neuron_lac_L_in (neuronas): ', num2str(lactateFluxNeuron)]);

% Verificar el flujo de lactato exportado por los astrocitos
lactateFluxAstro = solution.x(find(strcmp(mergedModel.rxns, 'EX_lac_L[e]')));
disp(['Flujo de lactato en la reacción EX_lac_L[e] (astrocitos): ', num2str(lactateFluxAstro)]);

% Verificar flujos nuevamente
gabaFluxAstro = solution.x(astroGabaIndex);
disp(['Flujo de GABA en la reacción Astro_gg4abut_in (astrocitos): ', num2str(gabaFluxAstro)]);

% Verificar el flujo de GABA liberado por las neuronas
gabaFluxNeuron = solution.x(find(strcmp(mergedModel.rxns, 'EX_gg4abut[e]')));
disp(['Flujo de GABA en la reacción EX_gg4abut[e] (neuronas): ', num2str(gabaFluxNeuron)]);

% Verificar el flujo de glutamato captado por los astrocitos
glutamateFluxAstro = solution.x(find(strcmp(mergedModel.rxns, 'Astro_glu_L_in')));
disp(['Flujo de glutamato en la reacción Astro_glu_L_in (astrocitos): ', num2str(glutamateFluxAstro)]);

% Verificar el flujo de glutamato exportado por las neuronas
glutamateFluxNeuron = solution.x(find(strcmp(mergedModel.rxns, 'EX_glu_L[e]')));
disp(['Flujo de glutamato en la reacción EX_glu_L[e] (neuronas): ', num2str(glutamateFluxNeuron)]);

% Verificar el flujo de alanina captada por los astrocitos
alanineFluxAstro = solution.x(find(strcmp(mergedModel.rxns, 'Astro_ala_D_in')));
disp(['Flujo de alanina D en la reacción Astro_ala_D_in (astrocitos): ', num2str(alanineFluxAstro)]);

solution = optimizeCbModel(mergedModel);

% Verificar los flujos nuevamente
aspFluxAstro = solution.x(find(strcmp(mergedModel.rxns, 'Neuron_asp_D_in')));
disp(['Flujo de aspartato en la reacción Neuron_asp_D_in (astrocitos): ', num2str(aspFluxAstro)]);
aspFluxNeuron = solution.x(aspExportIndex);
disp(['Flujo de aspartato en la reacción EX_asp_D[e] (neuronas): ', num2str(aspFluxNeuron)]);



%%%%_________


%%%%%%%%%%%_____________________________________



% Ajustar los límites de las reacciones de intercambio
% Por ejemplo, permitir que el lactato fluya del astrocito al compartimento extracelular
mergedModel = changeRxnBounds(mergedModel, 'EX_lac_L[e]', -10, 'l');  % Limitar la exportación de lactato
mergedModel = changeRxnBounds(mergedModel, 'Neuron_lac_L_in', 10, 'u');  % Permitir la entrada de lactato en neuronas

% Puedes ajustar otras reacciones de intercambio clave de manera similar

%_____________________________
% Añadir reacciones de intercambio para lactato, glutamato y GABA entre astrocito y neurona
% Lactato (producido por astrocitos, consumido por neuronas)
