import time
import os
import multiprocessing
import pickle
import requests


SHOW_PROGRESS_INTERVAL = 10
TEMP_FILE = f'static/ensembl_temp'
CHUNK_SIZE = 500
START = 15000
TOTAL_LENGTH = 18000
WAIT_INTERVAL = 0.2 # should limit to ~5 requests/second


from ccd.sequence_classes import HomologSequence

#Adding a counter to deal with anxiety. Adapted from
#https://stackoverflow.com/questions/2080660/python-multiprocessing-and-a-shared-counter

class Counter(object):
    
    def __init__(self, start=0):
        self.val = multiprocessing.Value('i', start)

    def next_(self, n=1):
        with self.val.get_lock():
            self.val.value += n
        return self.val.value

    @property
    def value(self):
        return self.val.value
    
def init_counter(counter):
    global global_counter
    global_counter = counter

def parse_ensembl_alignment(alignment, uniprot_id = None):
    sequences = []
    for hit in alignment['data'][0]['homologies']:
        species = hit['target']['species'] 
        protein_seq = hit['target']['align_seq']#.replace('-', '') -> already aligned?
        ensembl_id_ = hit['target']['id']
        #return namedtuple?
        homology_type = hit['type']
        sequences.append(HomologSequence(ensembl_id_, protein_seq, species,
                                         homology_type=homology_type))
    return sequences

def retrieve_ensembl_alignment(ensembl_id):
    if ensembl_id is None:
        return
    url = f'https://rest.ensemblgenomes.org//homology/id/{ensembl_id}'
    headers = {"Content-Type" : "application/json",
               "User-agent": "ccd - ccd@gmail.com"}
    response = requests.get(url, headers=headers)
    time.sleep(WAIT_INTERVAL) #to avoid flooding ensembl
    if not response.ok:
        return ('HttpError', response)
    else:
        return ('Success', response)

def get_alignment_from_ensembl(uniprot_id, ensembl_id):
#     print(f'Called with {uniprot_id}, {ensembl_id}')
    #response == ()
    if ensembl_id is None:
        return ('IdIsNone', None)
    status_code, requests_response = retrieve_ensembl_alignment(ensembl_id)
    if status_code == 'Success':
        try:
            alignment = parse_ensembl_alignment(requests_response.json())
            return (status_code, alignment)
        except Exception as e:
            print(f'exception on {uniprot_id}: {type(e).__name__}')
            return ('PythonError', e, requests_response)
    elif status_code == 'HttpError':
        print(f'{uniprot_id} failed due to {requests_response.status_code}')
        return (status_code, requests_response)

def split_work(uniprot_id, ensembl_id, queue): #uniprot_id, ensembl_id, queue
    global global_counter
    local_counter = global_counter.next_()
    print(f'working on entry {local_counter}')
    res = (uniprot_id, get_alignment_from_ensembl(uniprot_id, ensembl_id)) #tuple
    if not local_counter % SHOW_PROGRESS_INTERVAL:
        print(f'Processed entry {local_counter}')
    queue.put(res)
    

def listener(queue):
    print('I have been initialized')
    '''listens for messages on the queue, writes to file. '''
    with open(TEMP_FILE, 'wb') as f: 
        print(f'I opened the file')
        while not queue.empty():
            print('I am doing shtuff')
            res = queue.get_nowait()
            if res == 'kill':
                break
            pickle.dump(res) #pickle appends by default
            f.flush()

#processed to 7500-ish. before blocking. Check and continue
if __name__ == '__main__':
    multiprocessing.freeze_support()
    ensembl_references_file = r'static/ensembl_references.pickle'
    ensembl_alignments_file = r'static/ensembl_alignments.pickle'
    missing_ensembl_file = r'static/missing_ensembl.pickle'
    errors_file = r'static/failed_ensembl.pickle'
    save_files = False
    try:
        with open(ensembl_references_file, 'rb') as f:
            ensembl_references = pickle.load(f)
            print('Input file Read')
    except OSError:
        raise OSError (f'cannot find a file containing pdb_crossreferences at {ensembl_references_file}; please re run the notebook')
    global_counter = Counter(start=START)
    results = {}
    terminate = False
    if TOTAL_LENGTH: #just in case we want to retrieve only a chunk
        total_length = TOTAL_LENGTH
    else:
        total_length = len(ensembl_references) 
    for start in range(START, total_length, CHUNK_SIZE):
        stop = start + CHUNK_SIZE
        if stop >= total_length:
            terminate = True
            stop = total_length + 1
        with multiprocessing.Pool(processes= os.cpu_count() +1,
                          initializer=init_counter,
                          initargs=(global_counter,)) as pool:
            m = multiprocessing.Manager()
            queue = m.JoinableQueue()
            reference_list = [(*i, queue) for i in list(ensembl_references.items())[start:stop]]
            pool.starmap(split_work, reference_list)
            pool.close()
            pool.join()
            results = []
            while not queue.empty():
                results.append(queue.get())
            with open(f'{TEMP_FILE}_{start}_{stop}.pickle', 'wb') as f:
                pickle.dump(results, f)
                time.sleep(30) #give ENSEMBL some time
        if terminate:
            break
    thawed = []
    with open(f'{TEMP_FILE}_{start}_{stop}.pickle','rb') as f:
        while 1:
            try:
                thawed.append(pickle.load(f))
            except EOFError:
                break
    print(thawed) 
    
            
# def blocking_io():
#     # File operations (such as logging) can block the
#     # event loop: run them in a thread pool.
#     with open(os.path.abspath(r'static/uniprot.xml'), 'rb') as f:
#         return f.read(1000)
# 
# def cpu_bound():
#     # CPU-bound operations will block the event loop:
#     # in general it is preferable to run them in a
#     # process pool.
#     return sum(i * i for i in range(10 ** 7))
# async def main():
#     loop = asyncio.get_event_loop()
#     ## Options:
# 
#     # 1. Run in the default loop's executor:
#     result = await loop.run_in_executor(
#         None, blocking_io)
#     print('default thread pool', result)
# 
#     # 2. Run in a custom thread pool:
#     with concurrent.futures.ThreadPoolExecutor() as pool:
#         result = await loop.run_in_executor(
#             pool, blocking_io)
#         print('custom thread pool', result)
# 
#     # 3. Run in a custom process pool:
#     with concurrent.futures.ProcessPoolExecutor() as pool:
#         result = await loop.run_in_executor(
#             pool, cpu_bound)
#         print('custom process pool', result)
    
#     loop = asyncio.get_event_loop()
#     loop.run_until_complete(main())