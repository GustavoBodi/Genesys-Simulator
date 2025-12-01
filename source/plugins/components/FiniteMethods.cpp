/*
 * File:   FiniteMethods.cpp
 * Author: codex
 */

#include "FiniteMethods.h"
#include <fstream>
#include <list>
#include "../../kernel/simulator/Model.h"
#include "../../kernel/simulator/SimulationControlAndResponse.h"

#ifdef PLUGINCONNECT_DYNAMIC
extern "C" StaticGetPluginInformation GetPluginInformation() {
	return &FiniteMethods::GetPluginInformation;
}
#endif

ModelDataDefinition* FiniteMethods::NewInstance(Model* model, std::string name) {
	return new FiniteMethods(model, name);
}

std::string FiniteMethods::convertEnumToStr(DiffusionModel model) {
	switch (static_cast<int>(model)) {
		case 0: return "NONE";
		case 1: return "CONSTANT";
		default: return "UNKNOWN";
	}
}

std::string FiniteMethods::convertEnumToStr(Discretization discretization) {
	switch (static_cast<int>(discretization)) {
		case 0: return "FORWARD_EULER";
		case 1: return "CRANK_NICOLSON";
		default: return "UNKNOWN";
	}
}

std::string FiniteMethods::convertEnumToStr(Solver solver) {
	switch (static_cast<int>(solver)) {
		case 0: return "EXPLICIT";
		case 1: return "GAUSS_SEIDEL";
		default: return "UNKNOWN";
	}
}

FiniteMethods::FiniteMethods(Model* model, std::string name) : ModelComponent(model, Util::TypeOf<FiniteMethods>(), name) {
	SimulationControlGenericClass<Variable*, Model*, Variable>* propFieldVariable = new SimulationControlGenericClass<Variable*, Model*, Variable>(
									_parentModel,
									std::bind(&FiniteMethods::getFieldVariable, this), std::bind(&FiniteMethods::setFieldVariable, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "FieldVariable", "");
	SimulationControlGenericClass<Variable*, Model*, Variable>* propTimeVariable = new SimulationControlGenericClass<Variable*, Model*, Variable>(
									_parentModel,
									std::bind(&FiniteMethods::getTimeVariable, this), std::bind(&FiniteMethods::setTimeVariable, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "TimeVariable", "");
	SimulationControlGenericEnum<FiniteMethods::DiffusionModel, FiniteMethods>* propDiffusionModel = new SimulationControlGenericEnum<FiniteMethods::DiffusionModel, FiniteMethods>(
									std::bind(&FiniteMethods::getDiffusionModel, this), std::bind(&FiniteMethods::setDiffusionModel, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "DiffusionModel", "");
	SimulationControlGenericEnum<FiniteMethods::Discretization, FiniteMethods>* propDiscretization = new SimulationControlGenericEnum<FiniteMethods::Discretization, FiniteMethods>(
									std::bind(&FiniteMethods::getDiscretization, this), std::bind(&FiniteMethods::setDiscretization, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "Discretization", "");
	SimulationControlGenericEnum<FiniteMethods::Solver, FiniteMethods>* propSolver = new SimulationControlGenericEnum<FiniteMethods::Solver, FiniteMethods>(
									std::bind(&FiniteMethods::getSolver, this), std::bind(&FiniteMethods::setSolver, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "Solver", "");
	SimulationControlGeneric<double>* propDiffusionCoefficient = new SimulationControlGeneric<double>(
									std::bind(&FiniteMethods::getDiffusionCoefficient, this), std::bind(&FiniteMethods::setDiffusionCoefficient, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "DiffusionCoefficient", "");
	SimulationControlGeneric<double>* propSpaceStep = new SimulationControlGeneric<double>(
									std::bind(&FiniteMethods::getSpaceStep, this), std::bind(&FiniteMethods::setSpaceStep, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "SpaceStep", "");
	SimulationControlGeneric<double>* propTimeStep = new SimulationControlGeneric<double>(
									std::bind(&FiniteMethods::getTimeStep, this), std::bind(&FiniteMethods::setTimeStep, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "TimeStep", "");
	SimulationControlGeneric<std::string>* propFilename = new SimulationControlGeneric<std::string>(
									std::bind(&FiniteMethods::getFilename, this), std::bind(&FiniteMethods::setFilename, this, std::placeholders::_1),
									Util::TypeOf<FiniteMethods>(), getName(), "FileName", "");

	_parentModel->getControls()->insert(propFieldVariable);
	_parentModel->getControls()->insert(propTimeVariable);
	_parentModel->getControls()->insert(propDiffusionModel);
	_parentModel->getControls()->insert(propDiscretization);
	_parentModel->getControls()->insert(propSolver);
	_parentModel->getControls()->insert(propDiffusionCoefficient);
	_parentModel->getControls()->insert(propSpaceStep);
	_parentModel->getControls()->insert(propTimeStep);
	_parentModel->getControls()->insert(propFilename);

	_addProperty(propFieldVariable);
	_addProperty(propTimeVariable);
	_addProperty(propDiffusionModel);
	_addProperty(propDiscretization);
	_addProperty(propSolver);
	_addProperty(propDiffusionCoefficient);
	_addProperty(propSpaceStep);
	_addProperty(propTimeStep);
	_addProperty(propFilename);
}

void FiniteMethods::setFieldVariable(Variable* field) {
	_fieldVariable = field;
}

Variable* FiniteMethods::getFieldVariable() const {
	return _fieldVariable;
}

void FiniteMethods::setTimeVariable(Variable* timeVariable) {
	_timeVariable = timeVariable;
}

Variable* FiniteMethods::getTimeVariable() const {
	return _timeVariable;
}

void FiniteMethods::setDiffusionModel(DiffusionModel model) {
	_diffusionModel = model;
}

FiniteMethods::DiffusionModel FiniteMethods::getDiffusionModel() const {
	return _diffusionModel;
}

void FiniteMethods::setDiscretization(Discretization discretization) {
	_discretization = discretization;
}

FiniteMethods::Discretization FiniteMethods::getDiscretization() const {
	return _discretization;
}

void FiniteMethods::setSolver(Solver solver) {
	_solver = solver;
}

FiniteMethods::Solver FiniteMethods::getSolver() const {
	return _solver;
}

void FiniteMethods::setDiffusionCoefficient(double diffusionCoefficient) {
	_diffusionCoefficient = diffusionCoefficient;
}

double FiniteMethods::getDiffusionCoefficient() const {
	return _diffusionCoefficient;
}

void FiniteMethods::setSpaceStep(double spaceStep) {
	_spaceStep = spaceStep;
}

double FiniteMethods::getSpaceStep() const {
	return _spaceStep;
}

void FiniteMethods::setTimeStep(double timeStep) {
	_timeStep = timeStep;
}

double FiniteMethods::getTimeStep() const {
	return _timeStep;
}

void FiniteMethods::setFilename(std::string filename) {
	_filename = filename;
}

std::string FiniteMethods::getFilename() const {
	return _filename;
}

std::string FiniteMethods::show() {
	return ModelComponent::show() + "";
}

ModelComponent* FiniteMethods::LoadInstance(Model* model, PersistenceRecord *fields) {
	FiniteMethods* newComponent = new FiniteMethods(model);
	try {
		newComponent->_loadInstance(fields);
	} catch (const std::exception& e) {
	}
	return newComponent;
}

bool FiniteMethods::_doStep() {
	if (_fieldVariable == nullptr || _timeVariable == nullptr) {
		return false;
	}
	double now = _parentModel->getSimulation()->getSimulatedTime();
	double currentTime = _timeVariable->getValue();
	if (currentTime + _timeStep > now + 1e-15) {
		return false;
	}
	std::list<unsigned int>* dims = _fieldVariable->getDimensionSizes();
	if (dims->empty()) {
		return false;
	}
	unsigned int points = dims->front();
	if (points < 3) {
		return false;
	}
	std::vector<double> state(points, 0.0);
	for (unsigned int i = 0; i < points; ++i) {
		state[i] = _fieldVariable->getValue(std::to_string(i));
	}
	std::vector<double> next = state;
	_applyDiffusion(state, next);
	for (unsigned int i = 0; i < points; ++i) {
		_fieldVariable->setValue(next[i], std::to_string(i));
	}
	_timeVariable->setValue(currentTime + _timeStep);
	return true;
}

void FiniteMethods::_applyDiffusion(const std::vector<double>& current, std::vector<double>& next) {
	if (_diffusionModel == DiffusionModel::NONE || current.size() < 3) {
		next = current;
		return;
	}
	double alpha = (_spaceStep > 0.0) ? (_diffusionCoefficient * _timeStep) / (_spaceStep * _spaceStep) : 0.0;
	if (_solver == Solver::EXPLICIT) {
		next = current;
		if (_discretization == Discretization::CRANK_NICOLSON) {
			std::vector<double> predictor = current;
			for (std::size_t i = 1; i + 1 < current.size(); ++i) {
				double laplace = current[i - 1] - 2 * current[i] + current[i + 1];
				predictor[i] = current[i] + alpha * laplace;
			}
			for (std::size_t i = 1; i + 1 < current.size(); ++i) {
				double laplaceNow = current[i - 1] - 2 * current[i] + current[i + 1];
				double laplacePred = predictor[i - 1] - 2 * predictor[i] + predictor[i + 1];
				next[i] = current[i] + 0.5 * alpha * (laplaceNow + laplacePred);
			}
		} else {
			for (std::size_t i = 1; i + 1 < current.size(); ++i) {
				double laplace = current[i - 1] - 2 * current[i] + current[i + 1];
				next[i] = current[i] + alpha * laplace;
			}
		}
	} else {
		next = current;
		for (std::size_t i = 1; i + 1 < current.size(); ++i) {
			double left = next[i - 1];
			double right = current[i + 1];
			double laplace = left - 2 * next[i] + right;
			double candidate = next[i] + alpha * laplace;
			if (_discretization == Discretization::CRANK_NICOLSON) {
				double laplacePred = left - 2 * candidate + right;
				candidate = next[i] + 0.5 * alpha * (laplace + laplacePred);
			}
			next[i] = candidate;
		}
	}
	if (!next.empty()) {
		next.front() = current.front();
		next.back() = current.back();
	}
}

void FiniteMethods::_onDispatchEvent(Entity* entity, unsigned int inputPortNumber) {
	std::ofstream savefile;
	if (_filename != "") {
		savefile.open(_filename, std::ofstream::app);
	}
	while (_doStep()) {
		if (_fieldVariable != nullptr) {
			std::string message = "t=" + std::to_string(_timeVariable->getValue());
			for (unsigned int i = 0; i < _fieldVariable->getDimensionSizes()->front(); ++i) {
				message += " ," + _fieldVariable->getName() + "[" + std::to_string(i) + "]=" + std::to_string(_fieldVariable->getValue(std::to_string(i)));
			}
			traceSimulation(this, message, TraceManager::Level::L8_detailed);
			if (savefile.is_open()) {
				message = std::to_string(_timeVariable->getValue());
				for (unsigned int i = 0; i < _fieldVariable->getDimensionSizes()->front(); ++i) {
					message += "\t" + std::to_string(_fieldVariable->getValue(std::to_string(i)));
				}
				savefile << message << std::endl;
			}
		}
	}
	if (savefile.is_open()) {
		savefile.close();
	}
	_parentModel->sendEntityToComponent(entity, getConnectionManager()->getFrontConnection());
}

bool FiniteMethods::_loadInstance(PersistenceRecord *fields) {
	bool res = ModelComponent::_loadInstance(fields);
	if (res) {
		_diffusionModel = static_cast<DiffusionModel>(fields->loadField("diffusionModel", static_cast<int>(DEFAULT.diffusionModel)));
		_discretization = static_cast<Discretization>(fields->loadField("discretization", static_cast<int>(DEFAULT.discretization)));
		_solver = static_cast<Solver>(fields->loadField("solver", static_cast<int>(DEFAULT.solver)));
		_diffusionCoefficient = fields->loadField("diffusionCoefficient", DEFAULT.diffusionCoefficient);
		_spaceStep = fields->loadField("spaceStep", DEFAULT.spaceStep);
		_timeStep = fields->loadField("timeStep", DEFAULT.timeStep);
		_filename = fields->loadField("filename", DEFAULT.filename);
	}
	return res;
}

void FiniteMethods::_saveInstance(PersistenceRecord *fields, bool saveDefaultValues) {
	ModelComponent::_saveInstance(fields, saveDefaultValues);
	fields->saveField("diffusionModel", static_cast<int>(_diffusionModel), static_cast<int>(DEFAULT.diffusionModel), saveDefaultValues);
	fields->saveField("discretization", static_cast<int>(_discretization), static_cast<int>(DEFAULT.discretization), saveDefaultValues);
	fields->saveField("solver", static_cast<int>(_solver), static_cast<int>(DEFAULT.solver), saveDefaultValues);
	fields->saveField("diffusionCoefficient", _diffusionCoefficient, DEFAULT.diffusionCoefficient, saveDefaultValues);
	fields->saveField("spaceStep", _spaceStep, DEFAULT.spaceStep, saveDefaultValues);
	fields->saveField("timeStep", _timeStep, DEFAULT.timeStep, saveDefaultValues);
	fields->saveField("filename", _filename, DEFAULT.filename, saveDefaultValues);
}

bool FiniteMethods::_check(std::string& errorMessage) {
	bool ok = true;
	if (_fieldVariable == nullptr) {
		errorMessage += "FieldVariable not set; ";
		ok = false;
	}
	if (_timeVariable == nullptr) {
		errorMessage += "TimeVariable not set; ";
		ok = false;
	}
	if (_spaceStep <= 0.0) {
		errorMessage += "SpaceStep must be positive; ";
		ok = false;
	}
	if (_timeStep <= 0.0) {
		errorMessage += "TimeStep must be positive; ";
		ok = false;
	}
	if (_diffusionModel == DiffusionModel::CONSTANT && _diffusionCoefficient < 0.0) {
		errorMessage += "DiffusionCoefficient must be non-negative; ";
		ok = false;
	}
	if (ok && _fieldVariable != nullptr && _fieldVariable->getDimensionSizes()->empty()) {
		errorMessage += "FieldVariable dimension not set; ";
		ok = false;
	}
	if (ok && _fieldVariable != nullptr && _fieldVariable->getDimensionSizes()->front() < 3) {
		errorMessage += "FieldVariable requires at least 3 points; ";
		ok = false;
	}
	if (ok && _filename != "") {
		try {
			std::ofstream savefile(_filename, std::ofstream::out);
			std::string header = _timeVariable != nullptr ? _timeVariable->getName() : "t";
			if (_fieldVariable != nullptr) {
				for (unsigned int i = 0; i < _fieldVariable->getDimensionSizes()->front(); ++i) {
					header += "\t" + _fieldVariable->getName() + "[" + std::to_string(i) + "]";
				}
			}
			savefile << header << std::endl;
			savefile.close();
		} catch (...) {
			errorMessage += "Could not open file for FiniteMethods output; ";
			ok = false;
		}
	}
	return ok;
}

PluginInformation* FiniteMethods::GetPluginInformation() {
	PluginInformation* info = new PluginInformation(Util::TypeOf<FiniteMethods>(), &FiniteMethods::LoadInstance, &FiniteMethods::NewInstance);
	info->setCategory("Continuous");
	return info;
}
