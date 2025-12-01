/*
 * File:   Smart_FiniteMethods.cpp
 * Author: codex
 */

#include "Smart_FiniteMethods.h"

#include "../../../../kernel/simulator/Simulator.h"
#include "../../../../plugins/components/Create.h"
#include "../../../../plugins/components/FiniteMethods.h"
#include "../../../../plugins/components/Dispose.h"
#include "../../../../plugins/data/Variable.h"
#include "../../../TraitsApp.h"

Smart_FiniteMethods::Smart_FiniteMethods() {
}

int Smart_FiniteMethods::main(int argc, char** argv) {
	Simulator* genesys = new Simulator();
	genesys->getTraceManager()->setTraceLevel(TraitsApp<GenesysApplication_if>::traceLevel);
	setDefaultTraceHandlers(genesys->getTraceManager());
	PluginManager* plugins = genesys->getPluginManager();
	plugins->autoInsertPlugins("autoloadplugins.txt");
	Model* model = genesys->getModelManager()->newModel();
	Create* create1 = plugins->newInstance<Create>(model);
	create1->setTimeBetweenCreationsExpression("0.25");
	Dispose* dispose1 = plugins->newInstance<Dispose>(model);

	Variable* field = plugins->newInstance<Variable>(model, "u");
	field->insertDimentionSize(5);
	field->setInitialValue(1.0, "0");
	field->setInitialValue(1.0, "4");
	field->setInitialValue(0.0, "1");
	field->setInitialValue(0.0, "2");
	field->setInitialValue(0.0, "3");

	Variable* time = plugins->newInstance<Variable>(model, "t");
	time->setInitialValue(0.0);

	FiniteMethods* finite = plugins->newInstance<FiniteMethods>(model);
	finite->setFieldVariable(field);
	finite->setTimeVariable(time);
	finite->setDiffusionModel(FiniteMethods::DiffusionModel::CONSTANT);
	finite->setDiffusionCoefficient(0.8);
	finite->setSpaceStep(0.25);
	finite->setTimeStep(0.05);
	finite->setDiscretization(FiniteMethods::Discretization::CRANK_NICOLSON);
	finite->setSolver(FiniteMethods::Solver::GAUSS_SEIDEL);
	finite->setFilename("./temp/Smart_FiniteMethods.outputdatafile.txt");

	create1->getConnectionManager()->insert(finite);
	finite->getConnectionManager()->insert(dispose1);

	model->getSimulation()->setReplicationLength(1.0);
	model->getSimulation()->setShowReportsAfterReplication(false);
	model->getSimulation()->setShowReportsAfterSimulation(false);
	model->save("./models/Smart_FiniteMethods.gen");
	model->getSimulation()->start();
	delete genesys;
	return 0;
}
