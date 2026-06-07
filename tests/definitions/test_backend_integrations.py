from __future__ import annotations

from molinspect import load
from molinspect.backends.interfaces import _parse_pisa_interfaces_xml
from molinspect.backends.interfaces import BiologicalInterfaceRecord, BiologicalInterfaceStore
from molinspect.backends.interactions import _definition_id_for_arpeggio_contacts
from molinspect.backends.pockets import PocketRecord, PocketStore, _parse_p2rank_predictions


def test_arpeggio_mapping_keeps_weak_polar_as_lower_signal_than_polar():
    assert _definition_id_for_arpeggio_contacts(("proximal", "weak_polar")) == "near"
    assert _definition_id_for_arpeggio_contacts(("vdw_clash", "polar")) == "polar_contact_candidate"
    assert _definition_id_for_arpeggio_contacts(("ionic", "vdw")) == "salt_bridge"


def test_p2rank_prediction_csv_maps_lining_residues_to_pocket_records(tmp_path):
    prediction_csv = tmp_path / "demo_predictions.csv"
    prediction_csv.write_text(
        "name,rank,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,"
        "residue_ids,surf_atom_ids\n"
        "pocket1,1,9.77,0.525,70,40,70.5274,83.4375,-11.5099,A_103 A_180,32 33\n"
    )

    records = _parse_p2rank_predictions(prediction_csv)

    assert len(records) == 1
    record = records[0]
    assert record.id == "pocket:p2rank:1"
    assert record.backend == "P2Rank"
    assert record.center_A == (70.5274, 83.4375, -11.5099)
    assert record.probability == 0.525
    assert record.residue_expressions == (
        "chain A and resid 103",
        "chain A and resid 180",
    )
    assert record.annotations["method"] == "p2rank_ligand_binding_site_prediction"


def test_pisa_interfaces_xml_maps_buried_residues_to_biological_interface(tmp_path):
    interfaces_xml = tmp_path / "interfaces.xml"
    interfaces_xml.write_text(
        """
<pisa_results>
  <pdb_entry>
    <interface>
      <id>1</id>
      <int_area>1427.503</int_area>
      <int_solv_en>-18.217</int_solv_en>
      <pvalue>0.095</pvalue>
      <stab_en>-28.586</stab_en>
      <h-bonds><n_bonds>20</n_bonds></h-bonds>
      <salt-bridges><n_bonds>4</n_bonds></salt-bridges>
      <molecule>
        <chain_id>A-2</chain_id>
        <class>Protein</class>
        <residues>
          <residue>
            <name>ALA</name>
            <seq_num>2</seq_num>
            <ins_code></ins_code>
            <bsa>12.5</bsa>
          </residue>
          <residue>
            <name>VAL</name>
            <seq_num>3</seq_num>
            <ins_code></ins_code>
            <bsa>0.0</bsa>
          </residue>
        </residues>
      </molecule>
      <molecule>
        <chain_id>B</chain_id>
        <class>Protein</class>
        <residues>
          <residue>
            <name>LYS</name>
            <seq_num>4</seq_num>
            <ins_code>A</ins_code>
            <bsa>8.1</bsa>
          </residue>
        </residues>
      </molecule>
    </interface>
  </pdb_entry>
</pisa_results>
"""
    )

    records = _parse_pisa_interfaces_xml(interfaces_xml)

    assert len(records) == 1
    record = records[0]
    assert record.id == "biological_interface:pisa:1"
    assert record.backend == "PISA"
    assert record.chains == ("A", "B")
    assert record.pisa_chain_ids == ("A-2", "B")
    assert record.interface_area_A2 == 1427.503
    assert record.solvation_energy_kcal_mol == -18.217
    assert record.stabilization_energy_kcal_mol == -28.586
    assert record.pvalue == 0.095
    assert record.residue_expressions == (
        "chain A and resid 2 and resname ALA",
        "chain B and resid 4A and resname LYS",
    )
    assert record.annotations["bond_counts"] == {
        "hydrogen_bonds": 20,
        "salt_bridges": 4,
    }
    assert record.annotations["method"] == "pisa_surface_and_assembly_interface_analysis"


def test_pisa_legacy_interface_xml_maps_buried_residues(tmp_path):
    interfaces_xml = tmp_path / "legacy_interface.xml"
    interfaces_xml.write_text(
        """
<INTERFACE>
  <INTERFACENO>1</INTERFACENO>
  <INTERFACESUMMARY>
    <STRUCTURE1>
      <SOLVENTAREA1><INTERFACEAREA>849.515</INTERFACEAREA></SOLVENTAREA1>
      <SOLVATIONENERGY1>
        <GAINCOMPLEXFORMATION>-6.57079</GAINCOMPLEXFORMATION>
        <PVALUE>0.185634</PVALUE>
      </SOLVATIONENERGY1>
    </STRUCTURE1>
    <STRUCTURE2>
      <SOLVENTAREA2><INTERFACEAREA>843.744</INTERFACEAREA></SOLVENTAREA2>
      <SOLVATIONENERGY2>
        <GAINCOMPLEXFORMATION>-4.95907</GAINCOMPLEXFORMATION>
        <PVALUE>0.44685</PVALUE>
      </SOLVATIONENERGY2>
    </STRUCTURE2>
  </INTERFACESUMMARY>
  <HYDROGENBONDS>
    <STRUCTURE><STRUCTURE1>D:HIS 116[ NE2]</STRUCTURE1></STRUCTURE>
    <STRUCTURE><STRUCTURE1>D:ARG  30[ NH1]</STRUCTURE1></STRUCTURE>
  </HYDROGENBONDS>
  <RESIDUES>
    <RESIDUE1>
      <RESIDUE>
        <STRUCTURE> D:VAL   1    </STRUCTURE>
        <BURIEDSURFACEAREA>0</BURIEDSURFACEAREA>
      </RESIDUE>
      <RESIDUE>
        <STRUCTURE> D:ARG  30    </STRUCTURE>
        <BURIEDSURFACEAREA>80.3998</BURIEDSURFACEAREA>
      </RESIDUE>
    </RESIDUE1>
    <RESIDUE2>
      <RESIDUE>
        <STRUCTURE> C:PHE 117    </STRUCTURE>
        <BURIEDSURFACEAREA>31.8091</BURIEDSURFACEAREA>
      </RESIDUE>
    </RESIDUE2>
  </RESIDUES>
</INTERFACE>
"""
    )

    records = _parse_pisa_interfaces_xml(interfaces_xml)

    assert len(records) == 1
    record = records[0]
    assert record.id == "biological_interface:pisa:1"
    assert record.chains == ("D", "C")
    assert record.interface_area_A2 == 846.63
    assert record.solvation_energy_kcal_mol == -11.53
    assert record.pvalue == 0.447
    assert record.residue_expressions == (
        "chain D and resid 30 and resname ARG",
        "chain C and resid 117 and resname PHE",
    )
    assert record.annotations["bond_counts"] == {"hydrogen_bonds": 2}
    assert record.annotations["xml_variant"] == "pdbe_pisa_legacy_interface"


def test_session_objects_exposes_pisa_biological_interface_landmarks(tiny_pdb, monkeypatch):
    record = BiologicalInterfaceRecord(
        id="biological_interface:pisa:1",
        name="pisa_biological_interface 1",
        backend="PISA",
        pisa_id="1",
        chains=("A",),
        pisa_chain_ids=("A",),
        interface_area_A2=125.0,
        solvation_energy_kcal_mol=-2.0,
        stabilization_energy_kcal_mol=-3.0,
        pvalue=0.1,
        residue_ids=("residue:A:1:ALA",),
        residue_expressions=("chain A and resid 1 and resname ALA",),
        annotations={
            "definition_id": "biological_interface",
            "definition_source": "backend",
            "reference_keys": ["pisa_docs"],
            "backend": "PISA",
            "method": "pisa_surface_and_assembly_interface_analysis",
            "participating_residues": ["residue:A:1:ALA"],
            "participating_residue_count": 1,
        },
    )

    def fake_detect_biological_interfaces(universe, source_files):
        return BiologicalInterfaceStore(backend="PISA", records=(record,))

    monkeypatch.setattr(
        "molinspect.landmarks.detect_biological_interfaces",
        fake_detect_biological_interfaces,
    )

    session = load(structure=tiny_pdb)
    objects = session.objects(type="biological_interface")
    located = session.locate("biological_interface:pisa:1")

    assert objects.count == 1
    assert objects.objects[0].id == "biological_interface:pisa:1"
    assert located.location.objects == ["residue:A:1:ALA"]
    assert located.location.biological_interface_ids == ["biological_interface:pisa:1"]


def test_ranked_backend_landmarks_sort_by_numeric_rank(tiny_pdb, monkeypatch):
    records = tuple(
        PocketRecord(
            id=f"pocket:p2rank:{rank}",
            name=f"pocket{rank}_p2rank",
            backend="P2Rank",
            rank=rank,
            center_A=None,
            score=None,
            probability=None,
            residue_ids=("residue:A:1:ALA",),
            residue_expressions=("chain A and resid 1 and resname ALA",),
            annotations={"backend": "P2Rank", "rank": rank},
        )
        for rank in (1, 10, 2)
    )

    def fake_detect_pockets(universe, source_files):
        return PocketStore(backend="P2Rank", records=records)

    monkeypatch.setattr("molinspect.landmarks.detect_pockets", fake_detect_pockets)

    session = load(structure=tiny_pdb)
    objects = session.objects(type="pocket", limit=10)

    assert [obj.id for obj in objects.objects] == [
        "pocket:p2rank:1",
        "pocket:p2rank:2",
        "pocket:p2rank:10",
    ]
