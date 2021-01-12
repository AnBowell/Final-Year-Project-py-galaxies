Merger graph HDF storage structure
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   
This page details the data format for the graph files produced by MEGA, which is a work in progress and subject to change. 

Any dataset of length N host , N sub , N gen and N graph is ordered such that the index is the corresponding ID (all
IDs start at 0). For example, to get data for graph 3 from an N graph length array ( arr ) you would simply
index it thusly: 

arr[3]

Any dataset of length N prog , N desc , N subprog , and N subdesc contain the data for progenitors and descendents
and are of ambiguous length. Due to this there are associated pointer (start_index) and length (nprog,
ndesc) arrays. This is also true when trying to extract individual generations from the N host and N sub arrays
using the generation start and length arrays. 

Data can be extracted from these using:

arr[start: start + length]


****
Root
****

+-------------------------------------------------+
| **Groups**                                      |
+=========+=======================================+
| **Key** | **Description**                       |
+---------+---------------------------------------+
| Header  | Contains Metadata for the simulation. |
+---------+---------------------------------------+
| N       | | The groups containing the data for  |
|         | | each individual graph.              |
|         | | N runs from 0 - N graph -1.         |
+---------+---------------------------------------+

+--------------------------------------------------------------------------------------------------+
| **Datasets**                                                                                     |
+=====================+===============================================+===========+================+
| **Key**             | **Description**                               | **Units** | | **Type**/    |
|                     |                                               |           | | **Shape**    |
+---------------------+-----------------------------------------------+-----------+----------------+
| graph_lengths       | The length of each graph.                     | None      | | int[]/       |
|                     |                                               |           | | (N graph , ) |
+---------------------+-----------------------------------------------+-----------+----------------+
| root_nparts         | | The number of particles in the most         | None      | | int[]/       |
|                     | | massive host halo in the root generation    |           | | (N graph , ) |
|                     | | of each graph (lowest redshift generation). |           |                |
+---------------------+-----------------------------------------------+-----------+----------------+
| nhalos_in_graph     | The number of halos in each graph.            | None      | | int[]/       |
|                     |                                               |           | | (N graph , ) |
+---------------------+-----------------------------------------------+-----------+----------------+
| sub_graph_lengths   | The length of all subhalo graphs.             | None      | | int[]/       |
|                     |                                               |           | | (N graph , ) |
+---------------------+-----------------------------------------------+-----------+----------------+
| sub_nhalos_in_graph | The number of subhalos in each graph.         | None      | | int[]/       |
|                     |                                               |           | | (N graph , ) |
+---------------------+-----------------------------------------------+-----------+----------------+

******
Header
******

+-----------------------------------------------------------+
| **Attributes**                                            |
+===========+===============================+=======+=======+
| Key       | Description                   | Units | Type  |
+-----------+-------------------------------+-------+-------+
| past_mass | The dark matter particle mass | M_sun | float |
+-----------+-------------------------------+-------+-------+



**********
Graphs (N)
**********

+----------------------------------------------------------------------------------------+
| **Attributes**                                                                         |
+---------------------+-------------------------------------------+-----------+----------+
| **Key**             | **Description**                           | **Units** | **Type** |
+---------------------+-------------------------------------------+-----------+----------+
| length              | | The length of this graph, i.e.          | None      | int      |
|                     | | the number of generations               |           |          |
+---------------------+-------------------------------------------+-----------+----------+
| root_mass           | | The number of particles in the          | None      | int      |
|                     | | most massive halo in the root           |           |          |
|                     | | generation (lowest redshift generation) |           |          |
+---------------------+-------------------------------------------+-----------+----------+
| nhalos_in_graph     | | The number of halos in this             | None      | int      |
|                     | | graph                                   |           |          |
+---------------------+-------------------------------------------+-----------+----------+
| sub_length          | | The length of this subhalo graph,       | None      | int      |
|                     | | i.e. the number of generations          |           |          |
|                     | | containing subhalos                     |           |          |
+---------------------+-------------------------------------------+-----------+----------+
| sub_root_mass       | | The number of particles in the most     | None      | int      |
|                     | | massive subhalo in the root generation  |           |          |
|                     | | (lowest redshift generation)            |           |          |
+---------------------+-------------------------------------------+-----------+----------+
| sub_nhalos_in_graph | | The number of subhalos in this          | None      | int      |
|                     | | graph                                   |           |          |
+---------------------+-------------------------------------------+-----------+----------+
	 
	 
+------------------------------------------------------------------------------------------------------------+
| **Datasets**                                                                                               |
+==============================+================================================+===========+================+
| **Key**                      | **Description**                                | **Units** | | **Type**/    |
|                              |                                                |           | | **Shape**    |
+------------------------------+------------------------------------------------+-----------+----------------+
| graph_halo_ids               | | The internal graph halo ID assigned to a     | None      | | int[] /      |
|                              | | halo. NOTE: These differ from the halo       |           | | (N host , )  |
|                              | | catalog. These run from 0 - N host with the  |           |                |
|                              | | index equal to the value.                    |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| halo_catalog_halo_ids        | | The halo catalog ID assigned to each         | None      | | int[] /      |
|                              | | halo.                                        |           | | (N host , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| snapshots                    | | The index of the snapshot in the snapshot    | None      | | int[] /      |
|                              | | text file dictated in the param file,        |           | | (N host , )  |
|                              | | for each halo.                               |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| redshifts                    | | The redshift for each halo in the            | None      | | float[] /    |
|                              | | graph.                                       |           | | (N host , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| generation_id                | | The ID (or number) associated with each      | None      | | int[] /      |
|                              | | generation. Counting starts from             |           | | (N gen, )    |
|                              | | the earliest snapshot.                       |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| nparts                       | | The number of dark matter particles          | None      | | int[] /      |
|                              | | in each halo.                                |           | | (N host , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| mean_pos                     | | The mean position of the particles           | cMpc      | | float[] /    |
|                              | | in the halo.                                 |           | | (N host , 3) |
+------------------------------+------------------------------------------------+-----------+----------------+
| generation_start_index       | | The starting index (pointer) for each        | None      | | int[] /      |
|                              | | host halo generation.                        |           | | (N gen , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| generation_length            | | The number of halos in each                  | None      | | int[] /      |
|                              | | generation.                                  |           | | (N gen , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| nprog                        | | The number of progenitors for each           | None      | | int[] /      |
|                              | | halo.                                        |           | | (N host , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| ndesc                        | | The number of descendents for each           | None      | | int[] /      |
|                              | | halo.                                        |           | | (N host , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| prog_start_index             | | The starting index (pointer) for each halo’s | None      | | int[] /      |
|                              | | entries in all progenitor halo arrays (i.e.  |           | | (N host , )  |
|                              | | direct_prog_ids,                             |           |                |
|                              | | direct_prog_contribution, etc.).             |           |                |
|                              | | Entries containing 2**30 have no             |           |                |
|                              | | descendants.                                 |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| desc_start_index             | | The starting index (pointer) for each halo’s | None      | | int[] /      |
|                              | | entries in all descendant halo arrays(i.e.   |           | | (N host , )  |
|                              | | direct_desc_ids,                             |           |                |
|                              | | direct_desc_contribution, etc.).             |           |                |
|                              | | Entries containing 2**30 have no             |           |                |
|                              | | descendants.                                 |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| direct_prog_ids              | | The progenitor halo IDs, extracted using     | None      | | int[] /      |
|                              | | prog_start_index and nprog.                  |           | | (N prog, )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| direct_desc_ids              | | The descendent halo IDs, extracted using     | None      | | int[] /      |
|                              | | desc_start_index and ndesc.                  |           | | (N desc, )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| direct_prog_contribution     | | The number of dark matter particles          | None      | | int[] /      |
|                              | | contributed **by** each direct progenitor    |           | | (N prog , )  |
|                              | | **to** the halo.                             |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| direct_desc_contribution     | | The number of dark matter particles          | None      | | int[] /      |
|                              | | contributed **to** each direct descendent    |           | | (N desc, )   |
|                              | | **from** the halo.                           |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_graph_halo_ids           | | The internal graph subhalo ID assigned to a  | None      | | int[] /      |
|                              | | subhalo. NOTE: These differ from the halo    |           | | (N sub , )   |
|                              | | catalog. These run from 0 - N sut with the   |           |                |
|                              | | index equal to the value.                    |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| subhalo_catalog_halo_ids     | | The subhalo catalog ID assigned              | None      | | int[] /      |
|                              | | to each subhalo.                             |           | | (N sub , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_snapshots                | | The index of the snapshot in the snapshot    | None      | | int[] /      |
|                              | | text file dictated in the param file,        |           | | (N sub , )   |
|                              | | for each subhalo.                            |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_redshifts                | | The redshift for each subhalo in the         | None      | | float[] /    |
|                              | | graph.                                       |           | | (N sub ,)    |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_generation_id            | | The ID (or number) associated with each      | None      | | int[] /      |
|                              | | generation. Counting starts from the         |           | | (N gen , )   |
|                              | | earliest snapshot.                           |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_nparts                   | | The number of dark matter particles          | None      | | int[] /      |
|                              | | in each halo.                                |           | | (N sub , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_mean_pos                 | | The mean position of the particles           | None      | | float[] /    |
|                              | | in the subhalo.                              |           | | (N sub , 3)  |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_generation_start_index   | | The starting index (pointer) for each        | None      | | int[] /      |
|                              | | subhalo generation.                          |           | | (N gen , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_generation_length        | | The number of subhalos in each               | None      | | int[] /      |
|                              | | generation                                   |           | | (N gen , )   |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_nprog                    | | The number of progenitors for each           | None      | | int[] /      |
|                              | | subhalo.                                     |           | | (N subp , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_ndesc                    | | The number of descendents for                | None      | | int[] /      |
|                              | | each subhalo.                                |           | | (N subd , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_prog_start_index         | | The starting index (pointer) for each        | None      | | int[] /      |
|                              | | subhalo’s entries in all progenitor subhalo  |           | | (N subp , )  |
|                              | | arrays (i.e. sub_direct_prog_ids,            |           |                |
|                              | | sub_direct_prog_contribution, etc.).         |           |                |
|                              | | Entries containing 2**30 have no             |           |                |
|                              | | descendants.                                 |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_desc_start_index         | | The starting index (pointer) for each        | None      | | int[] /      |
|                              | | subhalo’s entries in all descendant subhalo  |           | | (N subd , )  |
|                              | | arrays (i.e. sub_direct_desc_ids,            |           |                |
|                              | | sub_direct_desc_contribution, etc.).         |           |                |
|                              | | Entries containing 2**30 have no             |           |                |
|                              | | descendants.                                 |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_direct_prog_ids          | | The progenitor subhalo IDs, extracted using  | None      | | int[] /      |
|                              | | sub_prog_start_index and sub_nprog .         |           | | (N subp , )  |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_direct_desc_ids          | | The descendent subhalo IDs, extracted        | None      | | int[] /      |
|                              | | using sub_desc_start_index and               |           | | (N subd , )  |
|                              | | sub_ndesc.                                   |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_direct_prog_contribution | | The number of dark matter particles          | None      | | int[] /      |
|                              | | contributed **by** each direct progenitor    |           | | (N subp , )  |
|                              | | **to** the subhalo.                          |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+
| sub_direct_desc_contribution | | The number of dark matter particles          | None      | | int[] /      |
|                              | | contributed to each direct descendent from   |           | | (N subd , )  |
|                              | | the subhalo.                                 |           |                |
+------------------------------+------------------------------------------------+-----------+----------------+


	 
	 
	 










