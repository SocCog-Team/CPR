function [exp_order, setup] = randomise_conditions()

% Random seed
rng("shuffle")

% Randomise setup allocation
setup = randi(3:4);
disp(['1st subject --> Setup: PSY' num2str(setup)]);

% Randomise conditions
conditions = {'Neutral', 'Cooperation', 'Competition'};
index = randperm(length(conditions));
exp_order = conditions(index);

disp('---');
disp(['Experimental Order: ']);
for i = 1:3
    disp([num2str(i) ') ' exp_order{i}]);
end

end

