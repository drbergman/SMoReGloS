function M = performEvents(M,in)

% performs the assigned tumor events

for ord_ind=1:length(in.active_ind)

    j = in.active_ind(ord_ind);

    switch in.events(ord_ind)
        case 1 % transition
            %%
            phase = M.tumor(j,M.I.phase);
            if M.chemo_pars.dna_check(phase) && rand() < M.chemo_pars.arrest_prob(phase)
                if M.flags.arrest_is_death % the cell dies
                    M.tumor(j,M.I.event) = 2; % mark the cell for apoptosis
                    in = stopFutureEvents(in,j,ord_ind);
                else % the cell becomes arrested
                    M.tumor(j,M.I.phase) = M.cycle.(M.cycle.phase_names(phase)+"a");
                    M.tumor(j,M.I.is_arrested) = true;
                    if M.pars.move_rate~=0 % the cell can undergo apoptosis from the arrested state
                        in = stopFutureMoveEvents(in,j,ord_ind);
                    end
                end
                M.tracked.chemo_arrest(M.i,phase) = M.tracked.chemo_arrest(M.i,phase)+1;
                continue;
            end

            if phase==M.cycle.m % then this cell proliferates

                [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumor(j,M.I.ind),M.tumor(j,M.I.subs),M); % neighbor indices
                if M.setup.ndims==3
                    if (nnz(M.L(n_ind)) + (26-length(n_ind)))<=M.pars.occmax % check how many M.pars.neighbors are occupied
                        weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                        ind = find(sum(weights)*rand() < cumsum(weights),1);
                        rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                        M.tumor(end+1,:) = 0;
                        M.tumor(end,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations
                        M.tumor(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                        M.L(n_ind(ind)) = M.val.tum; % set value at lattice site
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                        M.tumor(end,M.I.phase) = M.cycle.g1;

                        M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                    else
                        M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                    end
                else
                    if (nnz(M.L(n_ind)) + (8-length(n_ind)))<=M.pars.occmax % check how many M.pars.neighbors are occupied
                        weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                        ind = find(sum(weights)*rand() < cumsum(weights),1);
                        rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                        M.tumor(end+1,:) = 0;
                        M.tumor(end,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations
                        M.tumor(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                        M.L(n_ind(ind)) = M.val.tum; % set value at lattice site
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                        M.tumor(end,M.I.phase) = M.cycle.g1;

                        M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                    else
                        M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                    end

                end
            else % then it just moves along the transition path
                if M.tumor(j,M.I.is_arrested)==true
                    M.tracked.arrest_recovery(M.i,phase) = M.tracked.arrest_recovery(M.i,phase)+1;
                end
                M.tumor(j,M.I.phase) = M.cycle.advancer(M.tumor(j,M.I.phase));
                M.tumor(j,M.I.is_arrested) = false;
            end

        case 2 % spontaneous apoptosis
            %%
            M.L(M.tumor(j,M.I.ind)) = M.val.tum_apop;
            M.tracked.tum_apop(M.i) = M.tracked.tum_apop(M.i)+1;

        case 3 % random movement
            %%
            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumor(j,M.I.ind),M.tumor(j,M.I.subs),M); % neighbor indices
            if any(M.L(n_ind)==0)
                weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                ind = find(sum(weights)*rand() < cumsum(weights),1);
                rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                M.tumor(j,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations

                M.L(M.tumor(j,M.I.ind)) = 0;

                M.tumor(j,M.I.ind) = n_ind(ind); % store new array-specific locations
                M.L(n_ind(ind)) = M.val.tum; % set value at lattice site based on original cell (this should work instead of above commented out line)
            end

        case 4 % arrested cell skipped some event

        otherwise
            %%
            error('should not do nothing')

    end % end of switch
end % end of j for

end

function in = stopFutureEvents(in,j,ord_ind)
in.events(ord_ind + find(in.active_ind(ord_ind+1:end)==j)) = 4;
end

function in = stopFutureMoveEvents(in,j,ord_ind)
in.events(ord_ind + find(in.active_ind(ord_ind+1:end)==j & in.events(ord_ind+1:end)==3)) = 4;
end