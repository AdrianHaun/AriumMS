function IntResults = FinalizeIntegrationOutput(peaks,IntResults,EIC,times)  
% Performs Integration of found Peaks and sorts the output variable
                Areas=zeros(size(peaks,1),1);
                EIC = full(EIC);
                for n=1:size(peaks,1)
                    Areas(n,1)=trapz(EIC(peaks(n,2):peaks(n,3)));
                end
                if isempty(peaks)==false
                    IntResults{1,1}=peaks(:,4);
                    IntResults{2,1}=[Areas,peaks(:,2),peaks(:,3)];
                    IntResults{3,1}=[peaks(:,1),times(peaks(:,1))];
                    IntResults{8,1}= EIC;
                end
                
    end