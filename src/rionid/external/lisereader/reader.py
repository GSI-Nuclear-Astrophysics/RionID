import numpy as np
from rionid.external.barion.amedata import AMEData
from re import sub

class LISEreader:
    def __init__(self, filename):
        ame = AMEData()
        ame.init_ame_db
        self.ame_data = ame.ame_table
        self._read(filename)

    def _read(self, filename):
        with open(filename,encoding='latin1') as f:
            lines = f.readlines()

        for (i, line) in enumerate(lines):
            if '[Calculations]' in line:
                file_start = i+1
    
        # finds point in line where separator changes from space to comma
        self.centre_index = len(lines[file_start].split()) - 1

        # splits lines and adds to list of all data
        self.data = [line.replace('+', '').split()[0:self.centre_index] 
                     + line.replace('=', '').split()[self.centre_index].split(',')
                     for line in lines[file_start:]]

    def get_index(self, name):
        returned_queries = []
        for i, element in enumerate(self.data):
            if element[0] + '+' + element[1] == name:
                returned_queries.append(i)
            elif element[0] + element[1] + '+' == name:
                returned_queries.append(i)
                
        if len(returned_queries) == 1:
            return returned_queries[0]
        
        elif len(returned_queries) > 1:
            print('get_index() returned multiple queries, returning first')
            return returned_queries[0]
        
        else:
            print(element[0] + element[1] + '+')
            raise ValueError('get_index() returned nothing. Check formatting (e.g. \"80Kr35+\")')

    def get_info(self, name):
        try:
            index = self.get_index(name)
            from_lise = [list(map(int, self.data[index][1:self.centre_index])),
                         float(self.data[index][self.centre_index])]
    
            # Extract element and mass number
            element = sub('[^a-zA-Z]', '', name)
            aux = sub('[a-zA-Z]', '', name)
            aux2 = aux.index('+')
            mass_number = int(aux[:aux2])
    
            # Match to AME data
            from_ame_candidates = [[line[6], line[5], line[4], line[3]]
                                   for line in self.ame_data
                                   if (element == line[6] and mass_number == line[5])]
    
            if not from_ame_candidates:
                raise ValueError(
                    f"Nuclide {element}-{mass_number} not found in AME database.\n"
                    f"Please check if the AME data needs to be updated, or verify the formatting/spelling (e.g., 'Os' vs 'OS')."
                )
    
            from_ame = from_ame_candidates[0]
            return from_ame + from_lise
    
        except Exception as e:
            print(f"Error in get_info({name}): {e}")
            raise


    def get_lise_all(self):
        data = np.array(self.data)
        return np.transpose([data[:, 0], data[:, 5], data[:, 6]])

    def get_info_all(self):
        #return [self.get_info(line[0] + '+' + line[1]) for line in self.data]
        result = []
        for line in self.data:
            if len(line) >= 2:
                ion_name = line[0] + '+' + line[1]
                try:
                    info = self.get_info(ion_name)
                    result.append(info)
                except Exception as e:
                    print(f"Skipping {ion_name} due to error: {e}")
        return result    

    def get_info_specific(self,param_index_list):
        return_list = [[LISEreader.float_check(line[i]) for i in param_index_list]
                       for line in self.data]
        return return_list

    def get_yield_sorted(self):
        name_yield = [[str(line[0]) + '+' + str(line[1]), LISEreader.float_check(line[7])] for line in self.data]
        sorted_by_yield = sorted(name_yield, key=lambda x: x[1])[::-1]
        return sorted_by_yield
        
    @staticmethod
    def float_check(value): #returns float if not string
        if value.replace('.','').replace('e-','').replace('e+','').isdigit():
            return float(value)
        else:
            return value