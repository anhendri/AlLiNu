function [] = ALLINU_EXAM()
%% INITIALISATION DES MATRICES A TRAITER
    MAT_NAME = {'lund_a','bfw398b','bcsstk19','1138_bus','bcsstm27'};
%% CREATION DU FICHIER D'ENTREE ET DES MATRICES EN .dat
    FCND = fopen('CONDITION.dat','w');
    FIPT = fopen('INPUT.dat','w');
    fprintf(FIPT,'%d\n',numel(MAT_NAME));
    for INDX = 1:numel(MAT_NAME)
        NAME = MAT_NAME{INDX};
        fprintf(FIPT,'%s\n',NAME);
        MTRX = mmread([NAME '.mtx']);
        [POSX,POSY,VALR] = find(MTRX);
        FILE = fopen([NAME '.dat'],'w');
        fprintf(FILE,'%d\t%d\t%d\n',size(MTRX), length(POSX));
        fprintf(FILE,'%d\t%d\t%e\n',[POSX,POSY,VALR]');
        fclose(FILE);
%% ECRITURE DES TAILLES, DES CONDITIONNEMENTS ET DES NORMES DANS UN FICHIER
        fprintf(FCND,'%d\t',size(MTRX,1));
        MTRX = full(MTRX);
        fprintf(FCND,'%e\t',cond(MTRX,1),cond(MTRX,2),cond(MTRX,inf),       ...
                            norm(MTRX,1),norm(MTRX,2),norm(MTRX,inf));
        fprintf(FCND,'\n');
    end
    fclose(FCND);
    fclose(FIPT);
%% LANCEMENT DU PROGRAMME PRINCIPAL
    system('./ALLINU_EXAM.out');
%% CREATION DU FICHIER CONTENANT LES ANALYSES DES ERREURS
    FSOL = fopen('COMPARE.dat','w');
    for INDX = 1:numel(MAT_NAME)
        NAME = MAT_NAME{INDX};
%% RECUPERATION DU VECTEUR b ET DE LA MATRICE A
        VCTR = load(['rhs_' NAME '.dat']);
        MTRX = mmread([NAME '.mtx']);
        SLTN = [load(['result_' NAME '.dat']) MTRX\VCTR];
        FILE = fopen(['result_' NAME '.dat'],'w');
%% AJOUT DE LA SOLUTION MATLAB DANS LE FICHIER DES RESULTATS
        for POSX = 1:size(VCTR)
            for POSY = 1:3
                if(SLTN(POSX,POSY)>0)
                    fprintf(FILE,'+');
                end
                fprintf(FILE,'%0.10f\t',SLTN(POSX,POSY));
            end
            fprintf(FILE,'\n'); 
        end 
        fclose(FILE);
%% CALCUL DES ERREURS ABSOLUES ET RELATIVES
        ERRA = (SLTN(:,1:2) - [SLTN(:,3) SLTN(:,3)]);
        ERRA = abs(ERRA(SLTN(:,3)~=0,:));
        ERRR = (SLTN(:,1:2) - [SLTN(:,3) SLTN(:,3)])./[SLTN(:,3) SLTN(:,3)];
        ERRR = abs(ERRR(SLTN(:,3)~=0,:));
%% ECRITURE DANS LE FICHIER DES ERREURS ABSOLUES ET RELATIVES
        fprintf(FSOL,'%e\t',max(ERRA), norm(ERRA(:,1))/size(MTRX,1),       ...
                            norm(ERRA(:,2))/size(MTRX,1), mean(ERRA),      ...
                            max(ERRR), norm(ERRR(:,1))/size(MTRX,1),       ...
                            norm(ERRR(:,2))/size(MTRX,1), mean(ERRR));
        fprintf(FSOL,'\n');
    end
    fclose(FSOL);                                                             % FERMETURE DU FICHIER
end
