function kcal = calcola_kcal(eta, sesso, peso)
    if eta <= 10
        if peso <= 15
            kcal = 1000 + (peso-5)*300/10;
        elseif peso <= 25
            kcal = 1300 + (peso-15)*300/10;
        else
            kcal = 1600 + (peso-25)*400/10;
        end
    elseif eta <= 18
        if sesso == 'M'
            if peso <= 50
                kcal = 1500 + (peso-35)*300/15;
            elseif peso <= 60
                kcal = 1800 + (peso-50)*200/10;
            else
                kcal = 2000 + (peso-60)*500/15;
            end
        else
            if peso <= 45
                kcal = 1200 + (peso-30)*300/15;
            elseif peso <= 55
                kcal = 1500 + (peso-45)*300/10;
            else
                kcal = 1800 + (peso-55)*400/10;
            end
        end
    else % eta > 18
        if sesso == 'M'
            if peso <= 80
                kcal = 2200 + (peso-60)*300/20;
            elseif peso <= 100
                kcal =  2500 + (peso-80)*200/20;
            else
                kcal =  2700 + (peso-100)*300/20;
            end
        else
            if peso <= 55
                kcal = 2000 + (peso-45)*300/10;
            elseif peso <= 65
                kcal = 2300 + (peso-55)*200/10;
            else
                kcal = 2500 + (peso-65)*300/10;
            end
        end
    end
end
